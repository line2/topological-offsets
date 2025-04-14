#include "AmipsOffsetProblem.hpp"

#include <Eigen/Core>
#include <wmtk/components/topological_offsets/internal/utils/Quadrics.hpp>
#include <wmtk/function/utils/amips.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/simplex/link_single_dimension.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/simplex/utils/SimplexComparisons.hpp>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/orient.hpp>

namespace wmtk::operations {
AmipsOffsetProblem::AmipsOffsetProblem(
    const simplex::Simplex& s,
    const attribute::MeshAttributeHandle& embedding_pos,
    const attribute::MeshAttributeHandle& offset_pos,
    const SimpleBVH::BVH& bvh,
    const double offset_distance)
{
    // get optimal position
    {
        // computation of optimal point goes here
        const Mesh& m = offset_pos.mesh();
        assert(m.top_simplex_type() == PrimitiveType::Triangle);

        auto pos_offset_acc = m.create_const_accessor<double>(offset_pos);
        const Eigen::Vector3d p0 = pos_offset_acc.const_vector_attribute(s.tuple());

        // laplacian smoothing
        Eigen::Vector3d p_laplace = Eigen::Vector3d::Zero();
        const simplex::SimplexCollection neighs =
            simplex::link_single_dimension(m, s, PrimitiveType::Vertex);
        for (const simplex::Simplex& neigh : neighs) {
            p_laplace += pos_offset_acc.const_vector_attribute(neigh.tuple());
        }
        p_laplace /= neighs.size();

        // project with quadrics
        logger().info("Outdated and not maintained code");
        // wmtk::components::internal::utils::Quadrics q(bvh, offset_pos, s, offset_distance);
        // m_p_opt = q.solve(p_laplace);
    }

    // fill cells
    {
        const Mesh& m = embedding_pos.mesh();
        assert(m.top_simplex_type() == PrimitiveType::Tetrahedron);

        auto accessor = m.create_const_accessor(embedding_pos.as<double>());

        const simplex::Simplex s_embedding = offset_pos.mesh().map_to_parent(s);
        const auto neighs = wmtk::simplex::top_dimension_cofaces(m, s_embedding);

        m_cells.reserve(neighs.size());

        for (simplex::Simplex neigh : neighs) {
            if (!m.is_ccw(neigh.tuple())) {
                // switch any local id but NOT the vertex
                neigh = simplex::Simplex(
                    m,
                    neigh.primitive_type(),
                    m.switch_tuple(neigh.tuple(), PrimitiveType::Edge));
            }
            assert(m.is_ccw(neigh.tuple()));

            const auto vertices = simplex::faces_single_dimension(m, neigh, PrimitiveType::Vertex);
            assert(vertices.size() == 4);
            assert(simplex::utils::SimplexComparisons::equal(
                m,
                vertices.simplex_vector()[0],
                s_embedding));


            std::array<double, 12> cell;
            for (size_t i = 0; i < 4; ++i) {
                const simplex::Simplex& v = vertices.simplex_vector()[i];
                const auto p = accessor.const_vector_attribute(v.tuple());
                cell[3 * i + 0] = p[0];
                cell[3 * i + 1] = p[1];
                cell[3 * i + 2] = p[2];
            }
            m_cells.emplace_back(cell);
        }

        // check if input is valid
        {
            Eigen::Vector3d x0(0, 0, 0);
            Eigen::Vector3d x1(m_cells[0][0], m_cells[0][1], m_cells[0][2]);
            assert(is_step_valid(x0, x1));
        }

        m_inv_num_cells = 1. / m_cells.size();
    }
}

auto AmipsOffsetProblem::initial_value() const -> TVector
{
    Eigen::Vector3d tmp;
    for (int64_t d = 0; d < 3; ++d) {
        tmp(d) = m_cells[0][d]; // m_cells[0] is the point we want to optimize
    }

    return tmp;
}

double AmipsOffsetProblem::value(const TVector& x)
{
    assert(x.size() == 3);

    double amips_energy = 0;
    for (auto c : m_cells) {
        c[0] = x[0];
        c[1] = x[1];
        c[2] = x[2];
        amips_energy += wmtk::function::utils::Tet_AMIPS_energy(c);
    }

    double offset_energy = (x - m_p_opt).squaredNorm();

    return m_w_amips * m_inv_num_cells * amips_energy + m_w_offset * offset_energy;
}

void AmipsOffsetProblem::gradient(const TVector& x, TVector& gradv)
{
    constexpr int64_t size = 3;
    gradv.resize(3);
    gradv.setZero();
    Eigen::Vector3d tmp;

    Eigen::Vector3d amips_grad;
    amips_grad.setZero();
    for (auto c : m_cells) {
        c[0] = x[0];
        c[1] = x[1];
        c[2] = x[2];
        wmtk::function::utils::Tet_AMIPS_jacobian(c, tmp);
        amips_grad += tmp;
    }

    Eigen::Vector3d offset_grad = 2 * (x - m_p_opt);

    gradv = m_w_amips * m_inv_num_cells * amips_grad + m_w_offset * offset_grad;
}

void AmipsOffsetProblem::hessian(const TVector& x, Eigen::MatrixXd& hessian)
{
    hessian.resize(3, 3);
    hessian.setZero();
    Eigen::Matrix3d tmp;

    Eigen::Matrix3d amips_hess;
    amips_hess.setZero();
    for (auto c : m_cells) {
        c[0] = x[0];
        c[1] = x[1];
        c[2] = x[2];
        wmtk::function::utils::Tet_AMIPS_hessian(c, tmp);
        amips_hess += tmp;
    }

    Eigen::Matrix3d offset_hess = 2 * Eigen::Matrix3d::Identity();

    hessian = m_w_amips * m_inv_num_cells * amips_hess + m_w_offset * offset_hess;
}

void AmipsOffsetProblem::solution_changed(const TVector& new_x) {}

bool AmipsOffsetProblem::is_step_valid(const TVector& x0, const TVector& x1)
{
    const Eigen::Vector3d p0 = x1;

    for (const auto& c : m_cells) {
        const Eigen::Vector3d p1(c[3], c[4], c[5]);
        const Eigen::Vector3d p2(c[6], c[7], c[8]);
        const Eigen::Vector3d p3(c[9], c[10], c[11]);

        if (wmtk::utils::wmtk_orient3d(p3, p0, p1, p2) <= 0) {
            return false;
        }
    }

    return true;
}
} // namespace wmtk::operations