#include "Scheduler.hpp"

#include <wmtk/attribute/TypedAttributeHandle.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/simplex/k_ring.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/simplex/link.hpp>
#include <wmtk/simplex/utils/tuple_vector_to_homogeneous_simplex_vector.hpp>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/random_seed.hpp>

#include <polysolve/Utils.hpp>

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"
#endif
#include <tbb/parallel_for.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

// #include <tbb/task_arena.h>
#include <atomic>
#include <cstdlib>


#include <algorithm>
#include <random>

namespace wmtk {

Scheduler::Scheduler() = default;
Scheduler::~Scheduler() = default;

SchedulerStats Scheduler::run_operation_on_all(operations::Operation& op)
{
    SchedulerStats res;
    std::vector<simplex::Simplex> simplices;
    // op.reserve_enough_simplices();

    const auto type = op.primitive_type();
    {
        POLYSOLVE_SCOPED_STOPWATCH("Collecting primitives", res.collecting_time, logger());

        const auto tups = op.mesh().get_all(type);
        simplices =
            wmtk::simplex::utils::tuple_vector_to_homogeneous_simplex_vector(op.mesh(), tups, type);
    }

    logger().debug("Executing on {} simplices", simplices.size());
    std::vector<std::pair<int64_t, double>> order;

    {
        POLYSOLVE_SCOPED_STOPWATCH("Sorting", res.sorting_time, logger());
        if (op.use_random_priority()) {
            std::mt19937 gen(utils::get_random_seed());

            std::shuffle(simplices.begin(), simplices.end(), gen);
        } else {
            order.reserve(simplices.size());
            for (int64_t i = 0; i < simplices.size(); ++i) {
                order.emplace_back(i, op.priority(simplices[i]));
            }

            std::stable_sort(order.begin(), order.end(), [](const auto& s_a, const auto& s_b) {
                return s_a.second < s_b.second;
            });
        }
    }


    {
        POLYSOLVE_SCOPED_STOPWATCH("Executing operation", res.executing_time, logger());

        size_t total_simplices = simplices.size();

        if (op.use_random_priority()) {
            for (const auto& s : simplices) {
                log(res, total_simplices);
                auto mods = op(s);
                if (mods.empty())
                    res.fail();
                else
                    res.succeed();
            }
        } else {
            for (const auto& o : order) {
                log(res, total_simplices);
                auto mods = op(simplices[o.first]);
                if (mods.empty())
                    res.fail();
                else
                    res.succeed();
            }
        }
        if (m_update_frequency) {
            res.print_update_log(total_simplices);
        }
    }

    // logger().debug(
    //     "Ran {} ops, {} succeeded, {} failed",
    //     res.number_of_performed_operations(),
    //     res.number_of_successful_operations(),
    //     res.number_of_failed_operations());

    m_stats += res;

    return res;
}

SchedulerStats Scheduler::run_operation_on_all_parallel_prefilter(operations::Operation& op)
{
    SchedulerStats res;
    std::vector<simplex::Simplex> simplices;

    const auto type = op.primitive_type();
    {
        POLYSOLVE_SCOPED_STOPWATCH("Collecting primitives", res.collecting_time, logger());

        const auto tups = op.mesh().get_all(type);
        simplices =
            wmtk::simplex::utils::tuple_vector_to_homogeneous_simplex_vector(op.mesh(), tups, type);
    }

    logger().debug(
        "Executing on {} simplices (parallel prefilter)",
        simplices.size());
    std::vector<std::pair<int64_t, double>> order;

    {
        POLYSOLVE_SCOPED_STOPWATCH("Sorting", res.sorting_time, logger());
        if (op.use_random_priority()) {
            std::mt19937 gen(utils::get_random_seed());

            std::shuffle(simplices.begin(), simplices.end(), gen);
        } else {
            order.reserve(simplices.size());
            for (int64_t i = 0; i < simplices.size(); ++i) {
                order.emplace_back(i, op.priority(simplices[i]));
            }

            std::stable_sort(order.begin(), order.end(), [](const auto& s_a, const auto& s_b) {
                return s_a.second < s_b.second;
            });
        }
    }

    // Phase 1: parallel read-only prefilter (validity + before-invariants).
    // The mesh must not be mutated during this phase; each candidate is
    // re-validated at execution time below.
    // detail probes: measure each invariant separately (no behavior change)
    if (!m_probe_invariants.empty()) {
        std::vector<std::atomic<int64_t>> counts(m_probe_invariants.size());
        for (auto& c : counts) c = 0;
        tbb::parallel_for(tbb::blocked_range<int64_t>(0, (int64_t)simplices.size()), [&](const auto& r) {
            for (int64_t i = r.begin(); i < r.end(); ++i) {
                for (size_t k = 0; k < m_probe_invariants.size(); ++k) {
                    if (!m_probe_invariants[k]->before(simplices[i])) {
                        ++counts[k];
                    }
                }
            }
        });
        for (size_t k = 0; k < counts.size(); ++k) {
            m_probe_detail_counts[k] = counts[k].load();
        }
    }

    std::vector<char> keep(simplices.size(), 0);
    {
        POLYSOLVE_SCOPED_STOPWATCH("Parallel prefilter", res.prefilter_time, logger());

        const int64_t n = simplices.size();
        if (op.use_random_priority()) {
            tbb::parallel_for(tbb::blocked_range<int64_t>(0, n), [&](const auto& r) {
                for (int64_t i = r.begin(); i < r.end(); ++i) {
                    keep[i] = op.prefilter(simplices[i]) ? char(1) : char(0);
                }
            });
        } else {
            const bool gated = (m_probe_invariant != nullptr);
            const bool backoff = gated && m_backoff.has_value();
            const int64_t backoff_epoch = backoff ? m_backoff->epoch : 0;
            std::atomic<int64_t> probe_fails{0};
            std::atomic<int64_t> other_fails{0};
            std::atomic<int64_t> backoff_skips{0};
            if (backoff) {
                // created once per pass; shared read-only across threads
                auto stamp_acc =
                    op.mesh().create_const_accessor<int64_t>(m_backoff->fail_stamp_handle);
                auto dirty_acc =
                    op.mesh().create_const_accessor<int64_t>(m_backoff->dirty_epoch_handle);
                tbb::parallel_for(tbb::blocked_range<int64_t>(0, n), [&](const auto& r) {
                    for (int64_t i = r.begin(); i < r.end(); ++i) {
                        const int64_t idx = order[i].first;
                        if (gated && !m_probe_invariant->before(simplices[idx])) {
                            ++probe_fails;
                            keep[idx] = char(0);
                            continue;
                        }
                        const Tuple& t = simplices[idx].tuple();
                        const int64_t fs = stamp_acc.const_scalar_attribute(t);
                        if (fs > 0) {
                            const int64_t d0 = dirty_acc.const_scalar_attribute(t);
                            const int64_t d1 = dirty_acc.const_scalar_attribute(
                                op.mesh().switch_tuple(t, PrimitiveType::Vertex));
                            if (std::max(d0, d1) < fs) {
                                ++backoff_skips;
                                keep[idx] = char(0);
                                continue;
                            }
                        }
                        if (!op.prefilter(simplices[idx])) {
                            ++other_fails;
                            keep[idx] = char(0);
                            auto w_acc =
                                op.mesh().create_accessor<int64_t>(m_backoff->fail_stamp_handle);
                            w_acc.scalar_attribute(t) = backoff_epoch;
                            continue;
                        }
                        keep[idx] = char(1);
                    }
                });
            } else if (gated) {
                tbb::parallel_for(tbb::blocked_range<int64_t>(0, n), [&](const auto& r) {
                    for (int64_t i = r.begin(); i < r.end(); ++i) {
                        const int64_t idx = order[i].first;
                        if (!m_probe_invariant->before(simplices[idx])) {
                            ++probe_fails;
                            keep[idx] = char(0);
                        } else if (!op.prefilter(simplices[idx])) {
                            ++other_fails;
                            keep[idx] = char(0);
                        } else {
                            keep[idx] = char(1);
                        }
                    }
                });
            } else {
                tbb::parallel_for(tbb::blocked_range<int64_t>(0, n), [&](const auto& r) {
                    for (int64_t i = r.begin(); i < r.end(); ++i) {
                        const int64_t idx = order[i].first;
                        keep[idx] = op.prefilter(simplices[idx]) ? char(1) : char(0);
                    }
                });
            }
            m_probe_fails += probe_fails.load();
            m_probe_other_fails += other_fails.load();
            m_backoff_skips += backoff_skips.load();
        }
    }

    // Phase 2: serial execution of the survivors in (priority/random) order.
    {
        POLYSOLVE_SCOPED_STOPWATCH("Executing operation", res.executing_time, logger());

        size_t total_simplices = simplices.size();

        int64_t fails_since_last_success = 0;
        int64_t executed = 0;
        bool early_stopped = false;

        // Mark all vertices incident to a top-dim simplex containing the
        // operated edge (graph distance <= 1) dirty. This covers the
        // dependency set of all invariants used here:
        //  - energy/inversion/ROI invariants depend on the tets around the
        //    candidate edge, i.e. these vertices' positions;
        //  - the link condition of a candidate edge changes only when a
        //    tet containing one of its endpoints is modified, which implies
        //    the operated edge shared that tet, so the endpoint was marked;
        //  - moved/created vertices of the operation itself are endpoints of
        //    the operated edge, hence in the marked set.
        auto mark_neighborhood_dirty = [&](const simplex::Simplex& s) {
            if (!m_backoff.has_value()) return;
            m_backoff->epoch += 1;
            const int64_t epoch = m_backoff->epoch;
            auto dirty_acc = op.mesh().create_accessor<int64_t>(m_backoff->dirty_epoch_handle);
            for (const Tuple& ct : simplex::top_dimension_cofaces_tuples(op.mesh(), s)) {
                const simplex::Simplex cs(op.mesh(), op.mesh().top_simplex_type(), ct);
                const auto vs = simplex::faces_single_dimension_tuples(
                    op.mesh(),
                    cs,
                    PrimitiveType::Vertex);
                for (const Tuple& v : vs) {
                    dirty_acc.scalar_attribute(v) = epoch;
                }
            }
        };

        // Early-stop backoff: only kicks in for passes whose success rate has
        // collapsed (< 5%) and only after a long run of consecutive failures;
        // remaining candidates are retried by later passes.
        auto should_stop = [&]() {
            return m_early_stop > 0 && fails_since_last_success >= m_early_stop &&
                   executed >= 4 * m_early_stop &&
                   res.number_of_successful_operations() * 20 < executed;
        };

        if (op.use_random_priority()) {
            for (int64_t i = 0; i < (int64_t)simplices.size(); ++i) {
                if (keep[i] == char(0)) {
                    res.fail();
                    res.prefiltered();
                    continue;
                }
                if (early_stopped) {
                    res.fail();
                    continue;
                }
                log(res, total_simplices);
                auto mods = op(simplices[i]);
                ++executed;
                if (mods.empty()) {
                    res.fail();
                    ++fails_since_last_success;
                    if (should_stop()) {
                        early_stopped = true;
                    }
                } else {
                    res.succeed();
                    fails_since_last_success = 0;
                    // the operated simplex itself may no longer exist (e.g.
                    // after a split); mark the neighborhood of the resulting
                    // simplices instead
                    for (const auto& ms : mods) {
                        if (ms.primitive_type() == op.primitive_type()) {
                            mark_neighborhood_dirty(ms);
                        }
                    }
                }
            }
        } else {
            for (const auto& o : order) {
                if (keep[o.first] == char(0)) {
                    res.fail();
                    res.prefiltered();
                    continue;
                }
                if (early_stopped) {
                    res.fail();
                    continue;
                }
                log(res, total_simplices);
                auto mods = op(simplices[o.first]);
                ++executed;
                if (mods.empty()) {
                    res.fail();
                    ++fails_since_last_success;
                    if (should_stop()) {
                        early_stopped = true;
                    }
                } else {
                    res.succeed();
                    fails_since_last_success = 0;
                    for (const auto& ms : mods) {
                        if (ms.primitive_type() == op.primitive_type()) {
                            mark_neighborhood_dirty(ms);
                        }
                    }
                }
            }
        }
        if (m_update_frequency) {
            res.print_update_log(total_simplices);
        }
    }

    m_stats += res;

    return res;
}

SchedulerStats Scheduler::run_operation_on_all(
    operations::Operation& op,
    const TypedAttributeHandle<char>& flag_handle)
{
    std::vector<simplex::Simplex> simplices;
    const auto type = op.primitive_type();

    auto flag_accessor = op.mesh().create_accessor(flag_handle);
    auto tups = op.mesh().get_all(type);

    SchedulerStats res;
    int64_t success = -1;
    std::vector<std::pair<int64_t, double>> order;

    for (const auto& t : tups) {
        flag_accessor.scalar_attribute(t) = char(1);
    }

    do {
        SchedulerStats internal_stats;
        // op.reserve_enough_simplices();

        {
            POLYSOLVE_SCOPED_STOPWATCH(
                "Collecting primitives",
                internal_stats.collecting_time,
                logger());

            tups = op.mesh().get_all(type);
            const auto n_primitives = tups.size();
            tups.erase(
                std::remove_if(
                    tups.begin(),
                    tups.end(),
                    [&](const Tuple& t) { return flag_accessor.scalar_attribute(t) == char(0); }),
                tups.end());
            for (const auto& t : tups) {
                flag_accessor.scalar_attribute(t) = char(0);
            }

            logger().debug("Processing {}/{}", tups.size(), n_primitives);

            simplices = wmtk::simplex::utils::tuple_vector_to_homogeneous_simplex_vector(
                op.mesh(),
                tups,
                type);
        }

        {
            POLYSOLVE_SCOPED_STOPWATCH("Sorting", internal_stats.sorting_time, logger());
            if (op.use_random_priority()) {
                std::mt19937 gen(utils::get_random_seed());

                std::shuffle(simplices.begin(), simplices.end(), gen);
            } else {
                order.clear();
                order.reserve(simplices.size());
                for (int64_t i = 0; i < simplices.size(); ++i) {
                    order.emplace_back(i, op.priority(simplices[i]));
                }

                std::stable_sort(order.begin(), order.end(), [](const auto& s_a, const auto& s_b) {
                    return s_a.second < s_b.second;
                });
            }
        }

        {
            size_t total_simplices = order.size();
            POLYSOLVE_SCOPED_STOPWATCH(
                "Executing operation",
                internal_stats.executing_time,
                logger());

            if (op.use_random_priority()) {
                for (const auto& s : simplices) {
                    log(internal_stats, total_simplices);
                    auto mods = op(s);
                    // assert(wmtk::multimesh::utils::check_child_maps_valid(op.mesh()));
                    if (mods.empty()) {
                        res.fail();
                    } else {
                        res.succeed();
                    }
                }
            } else {
                for (const auto& o : order) {
                    log(internal_stats, total_simplices);
                    auto mods = op(simplices[o.first]);
                    // assert(wmtk::multimesh::utils::check_child_maps_valid(op.mesh()));
                    if (mods.empty()) {
                        internal_stats.fail();
                    } else {
                        internal_stats.succeed();
                    }
                }
            }
        }

        success = internal_stats.number_of_successful_operations();
        res += internal_stats;
        res.sub_stats.push_back(internal_stats);
        m_stats += internal_stats;
        m_stats.sub_stats.push_back(internal_stats);
    } while (success > 0);

    // // reset flag to 1, not necessaty
    // auto tups = op.mesh().get_all(type);
    // for (const auto& t : tups) {
    //     flag_accessor.scalar_attribute(t) = char(1);
    // }

    return res;
}

int64_t first_available_color(std::vector<int64_t>& used_neighbor_coloring)
{
    if (used_neighbor_coloring.size() == 0) return 0;

    std::sort(used_neighbor_coloring.begin(), used_neighbor_coloring.end());
    int64_t color = 0;

    for (const auto c : used_neighbor_coloring) {
        if (color < c) {
            return color;
        } else {
            color = c + 1;
        }
    }
    return color;
}

SchedulerStats Scheduler::run_operation_on_all_coloring(
    operations::Operation& op,
    const TypedAttributeHandle<int64_t>& color_handle,
    bool reuse_existing_colors)
{
    // this only works on vertex operations
    SchedulerStats res;
    std::vector<std::vector<simplex::Simplex>> colored_simplices;
    // op.reserve_enough_simplices();
    auto color_accessor = op.mesh().create_accessor<int64_t>(color_handle);

    const auto type = op.primitive_type();
    assert(type == PrimitiveType::Vertex);

    const auto tups = op.mesh().get_all(type);
    int64_t color_max = -1;
    if (!reuse_existing_colors) {
        POLYSOLVE_SCOPED_STOPWATCH("Collecting primitives", res.collecting_time, logger());

        tbb::parallel_for(
            tbb::blocked_range<int64_t>(0, tups.size()),
            [&](tbb::blocked_range<int64_t> r) {
                for (int64_t i = r.begin(); i < r.end(); ++i) {
                    color_accessor.scalar_attribute(tups[i]) = -1;
                }
            });

        // greedy coloring
        std::vector<int64_t> used_colors;
        for (const auto& v : tups) {
            // get used colors in neighbors
            used_colors.clear();
            for (const auto& v_one_ring :
                 simplex::link(op.mesh(), simplex::Simplex::vertex(op.mesh(), v), false)
                     .simplex_vector(PrimitiveType::Vertex)) {
                int64_t color = color_accessor.const_scalar_attribute(v_one_ring.tuple());
                if (color > -1) {
                    used_colors.push_back(color);
                }
            }
            int64_t c = first_available_color(used_colors);
            color_accessor.scalar_attribute(v) = c;
            color_max = std::max(color_max, c);

            // push into vectors
            if (c + 1 > colored_simplices.size()) {
                colored_simplices.push_back({simplex::Simplex::vertex(op.mesh(), v)});
            } else {
                colored_simplices[c].push_back(simplex::Simplex::vertex(op.mesh(), v));
            }
        }

        logger().info("Have {} colors among {} vertices", colored_simplices.size(), tups.size());

    } else {
        // bucket vertices by the caller-provided coloring
        for (const auto& v : tups) {
            const int64_t c = color_accessor.const_scalar_attribute(v);
            if (c < 0) {
                continue;
            }
            color_max = std::max(color_max, c);
            if (c + 1 > (int64_t)colored_simplices.size()) {
                colored_simplices.resize(c + 1);
            }
            colored_simplices[c].push_back(simplex::Simplex::vertex(op.mesh(), v));
        }
        logger().info(
            "Reusing existing coloring: {} colors among {} vertices",
            colored_simplices.size(),
            tups.size());
    }

    logger().debug("Executing on {} simplices", tups.size());

    {
        POLYSOLVE_SCOPED_STOPWATCH("Sorting", res.sorting_time, logger());

        // do sth or nothing
    }

    const char* serial_mode = std::getenv("WMTK_COLORING_SERIAL");
    bool force_serial = false;
    if (serial_mode != nullptr) {
        if (std::string(serial_mode) == "1") {
            force_serial = true;
        } else if (std::string(serial_mode) == "tet") {
            force_serial = (op.mesh().top_simplex_type() == PrimitiveType::Tetrahedron);
        } else if (std::string(serial_mode) == "tri") {
            force_serial = (op.mesh().top_simplex_type() == PrimitiveType::Triangle);
        }
    }

    // verify coloring validity (adjacent vertices must have different colors)
    if (std::getenv("WMTK_COLORING_CHECK") != nullptr) {
        int64_t conflicts = 0;
        for (const auto& v : tups) {
            const auto cv = color_accessor.const_scalar_attribute(v);
            for (const auto& v_one_ring :
                 simplex::link(op.mesh(), simplex::Simplex::vertex(op.mesh(), v), false)
                     .simplex_vector(PrimitiveType::Vertex)) {
                if (cv == color_accessor.const_scalar_attribute(v_one_ring.tuple())) {
                    ++conflicts;
                }
            }
        }
        if (conflicts > 0) {
            logger().error("COLORING CONFLICTS: {} adjacent same-color pairs", conflicts);
        } else {
            logger().info("coloring check passed");
        }
    }


    {
        POLYSOLVE_SCOPED_STOPWATCH("Executing operation", res.executing_time, logger());

        std::atomic_int suc_cnt = 0;
        std::atomic_int fail_cnt = 0;

        // same dirty-marking logic as the prefilter scheduler: on success,
        // mark all vertices sharing a top-dim simplex with the modified
        // simplex so that stamped failures in the neighborhood are
        // re-evaluated
        auto mark_dirty_coloring = [&](const simplex::Simplex& s) {
            if (!m_backoff.has_value()) return;
            m_backoff->epoch += 1;
            const int64_t epoch = m_backoff->epoch;
            auto dirty_acc = op.mesh().create_accessor<int64_t>(m_backoff->dirty_epoch_handle);
            for (const Tuple& ct : simplex::top_dimension_cofaces_tuples(op.mesh(), s)) {
                const simplex::Simplex cs(op.mesh(), op.mesh().top_simplex_type(), ct);
                const auto vs = simplex::faces_single_dimension_tuples(
                    op.mesh(),
                    cs,
                    PrimitiveType::Vertex);
                for (const Tuple& v : vs) {
                    dirty_acc.scalar_attribute(v) = epoch;
                }
            }
        };

        for (int64_t i = 0; i < colored_simplices.size(); ++i) {
            if (force_serial) {
                for (int64_t k = 0; k < (int64_t)colored_simplices[i].size(); ++k) {
                    auto mods = op(colored_simplices[i][k]);
                    if (mods.empty()) {
                        fail_cnt++;
                    } else {
                        suc_cnt++;
                        for (const auto& ms : mods) {
                            if (ms.primitive_type() == op.primitive_type()) {
                                mark_dirty_coloring(ms);
                            }
                        }
                    }
                }
                continue;
            }
            tbb::parallel_for(
                tbb::blocked_range<int64_t>(0, colored_simplices[i].size()),
                [&](tbb::blocked_range<int64_t> r) {
                    for (int64_t k = r.begin(); k < r.end(); ++k) {
                        auto mods = op(colored_simplices[i][k]);
                        if (mods.empty()) {
                            fail_cnt++;
                        } else {
                            suc_cnt++;
                            for (const auto& ms : mods) {
                                if (ms.primitive_type() == op.primitive_type()) {
                                    mark_dirty_coloring(ms);
                                }
                            }
                        }
                    }
                });

            // for (int64_t k = 0; k < colored_simplices[i].size(); ++k) {
            //     auto mods = op(colored_simplices[i][k]);
            //     if (mods.empty()) {
            //         fail_cnt++;
            //     } else {
            //         suc_cnt++;
            //     }
            // }
        }

        res.m_num_op_success = suc_cnt;
        res.m_num_op_fail = fail_cnt;
    }

    // logger().debug(
    //     "Ran {} ops, {} succeeded, {} failed",
    //     res.number_of_performed_operations(),
    //     res.number_of_successful_operations(),
    //     res.number_of_failed_operations());

    m_stats += res;

    return res;
}

void Scheduler::log(size_t total)
{
    log(m_stats, total);
}
void Scheduler::log(const SchedulerStats& stats, size_t total)
{
    if (m_update_frequency) {
        const size_t freq = m_update_frequency.value();
        size_t count = stats.number_of_performed_operations();
        if (count % freq == 0) {
            stats.print_update_log(total);
        }
        count++;
    }
}

void Scheduler::set_update_frequency(std::optional<size_t>&& freq)
{
    m_update_frequency = std::move(freq);
}


void SchedulerStats::print_update_log(size_t total, spdlog::level::level_enum level) const
{
    logger().log(
        level,
        "{} success + {} fail = {} out of {}",
        number_of_successful_operations(),
        number_of_failed_operations(),
        number_of_performed_operations(),
        total);
}

} // namespace wmtk
