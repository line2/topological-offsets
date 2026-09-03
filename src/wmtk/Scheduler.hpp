#pragma once

#include <spdlog/common.h>
#include <wmtk/operations/Operation.hpp>

namespace wmtk {

class SchedulerStats
{
public:
    /**
     * @brief Returns the number of successful operations performed by the scheduler.
     *
     * The value is reset to 0 when calling `run_operation_on_all`.
     */
    int64_t number_of_successful_operations() const { return m_num_op_success; }

    /**
     * @brief Returns the number of failed operations performed by the scheduler.
     *
     * The value is reset to 0 when calling `run_operation_on_all`.
     */
    int64_t number_of_failed_operations() const { return m_num_op_fail; }

    /**
     * @brief Returns the number of performed operations performed by the scheduler.
     *
     * The value is reset to 0 when calling `run_operation_on_all`.
     */
    int64_t number_of_performed_operations() const { return m_num_op_success + m_num_op_fail; }

    /**
     * @brief Operations rejected by the parallel prefilter without execution.
     *
     * These are also counted as failed operations for compatibility.
     */
    int64_t number_of_prefiltered_operations() const { return m_num_prefiltered; }

    inline double total_time() const
    {
        return collecting_time + sorting_time + prefilter_time + executing_time;
    }

    inline void succeed() { ++m_num_op_success; }
    inline void fail() { ++m_num_op_fail; }
    inline void prefiltered() { ++m_num_prefiltered; }

    inline void operator+=(const SchedulerStats& s)
    {
        m_num_op_success += s.m_num_op_success;
        m_num_op_fail += s.m_num_op_fail;
        m_num_prefiltered += s.m_num_prefiltered;

        collecting_time += s.collecting_time;
        sorting_time += s.sorting_time;
        prefilter_time += s.prefilter_time;
        executing_time += s.executing_time;
    }


    double collecting_time = 0;
    double sorting_time = 0;
    double prefilter_time = 0;
    double executing_time = 0;

    std::vector<SchedulerStats> sub_stats;

    double avg_sub_collecting_time() const
    {
        double res = 0;
        for (const auto& s : sub_stats) {
            res += s.collecting_time;
        }
        return res / sub_stats.size();
    }

    double avg_sub_sorting_time() const
    {
        double res = 0;
        for (const auto& s : sub_stats) {
            res += s.sorting_time;
        }
        return res / sub_stats.size();
    }

    double avg_sub_executing_time() const
    {
        double res = 0;
        for (const auto& s : sub_stats) {
            res += s.executing_time;
        }
        return res / sub_stats.size();
    }

    // private:
    int64_t m_num_op_success = 0;
    int64_t m_num_op_fail = 0;
    int64_t m_num_prefiltered = 0;

    void print_update_log(size_t total, spdlog::level::level_enum = spdlog::level::info) const;
};

class Scheduler
{
public:
    Scheduler();
    ~Scheduler();

    SchedulerStats run_operation_on_all(operations::Operation& op);
    SchedulerStats run_operation_on_all(
        operations::Operation& op,
        const TypedAttributeHandle<char>& flag_handle);

    /**
     * @brief Two-phase variant: parallel read-only prefilter + serial execution.
     *
     * Phase 1 evaluates validity and before-invariants for all candidate
     * simplices concurrently (read-only, mesh must not change). Phase 2
     * executes the survivors serially in priority order; each executed
     * operation re-validates against the current mesh state, so results are
     * consistent with the serial scheduler (some candidates may be rejected
     * a second time if earlier executions changed their neighborhood).
     */
    SchedulerStats run_operation_on_all_parallel_prefilter(operations::Operation& op);
    /**
     * @brief Executes a vertex operation in parallel over independent color
     * classes.
     *
     * @param reuse_existing_colors when true, the colors already stored in
     * color_handle are used as-is (vertices with color < 0 are skipped). This
     * allows callers to supply conflict colorings based on meshes other than
     * the operation's own mesh (e.g. coloring offset-surface vertices by
     * adjacency of their parents in the embedding mesh).
     */
    SchedulerStats run_operation_on_all_coloring(
        operations::Operation& op,
        const TypedAttributeHandle<int64_t>& color_handle,
        bool reuse_existing_colors = false);

    const SchedulerStats& stats() const { return m_stats; }

    void set_update_frequency(std::optional<size_t>&& freq = {});

    /**
     * @brief Early-stop backoff for the parallel-prefilter variant: stop the
     * serial execution phase once this many consecutively executed
     * operations fail (0 disables). Remaining candidates are counted as
     * failed and retried by later passes, mirroring prefilter rejections.
     */
    void set_early_stop_after_consecutive_failures(int64_t n) { m_early_stop = n; }

    /**
     * @brief Instrumentation: when set, the parallel prefilter evaluates this
     * invariant first and counts how many candidates fail it versus failing
     * the remaining before-invariants. Read via probe_fail_count()/
     * probe_other_fail_count() after the run.
     */
    void set_probe_invariant(std::shared_ptr<invariants::Invariant> inv)
    {
        m_probe_invariant = std::move(inv);
        m_probe_fails = 0;
        m_probe_other_fails = 0;
    }
    int64_t probe_fail_count() const { return m_probe_fails; }
    int64_t probe_other_fail_count() const { return m_probe_other_fails; }

    /**
     * @brief Instrumentation: evaluate each probe invariant on all
     * candidates (measurement only, does not change behavior) and count
     * individual failures. Read via probe_detail_counts().
     */
    void set_probe_invariants(std::vector<std::shared_ptr<invariants::Invariant>> invs)
    {
        m_probe_invariants = std::move(invs);
        m_probe_detail_counts.assign(m_probe_invariants.size(), 0);
    }
    const std::vector<int64_t>& probe_detail_counts() const { return m_probe_detail_counts; }

private:
    SchedulerStats m_stats;
    std::optional<size_t> m_update_frequency = {};
    int64_t m_early_stop = 0;
    std::shared_ptr<invariants::Invariant> m_probe_invariant = nullptr;
    int64_t m_probe_fails = 0;
    int64_t m_probe_other_fails = 0;
    std::vector<std::shared_ptr<invariants::Invariant>> m_probe_invariants;
    std::vector<int64_t> m_probe_detail_counts;

    void log(const size_t total);
    void log(const SchedulerStats& stats, const size_t total);
};

} // namespace wmtk
