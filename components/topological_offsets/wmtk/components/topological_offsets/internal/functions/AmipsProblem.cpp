#include "AmipsProblem.hpp"

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
AmipsProblem::AmipsProblem(
    const simplex::Simplex& s,
    const attribute::MeshAttributeHandle& pos_handle)
    : m_pt(pos_handle.mesh().top_simplex_type())
{
    // fill cells
    const Mesh& m = pos_handle.mesh();

    auto p_acc = m.create_const_accessor(pos_handle.as<double>());
    const auto neighs = wmtk::simplex::top_dimension_cofaces(m, s);

    if (m.top_simplex_type() == PrimitiveType::Tetrahedron) {
        m_cells_tet.reserve(neighs.size());

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
            assert(simplex::utils::SimplexComparisons::equal(m, vertices.simplex_vector()[0], s));


            std::array<double, 12> cell;
            for (size_t i = 0; i < 4; ++i) {
                const simplex::Simplex& v = vertices.simplex_vector()[i];
                const auto p = p_acc.const_vector_attribute(v.tuple());
                cell[3 * i + 0] = p[0];
                cell[3 * i + 1] = p[1];
                cell[3 * i + 2] = p[2];
            }
            m_cells_tet.emplace_back(cell);
        }

        // check if input is valid
        {
            Eigen::Vector3d x0(0, 0, 0);
            Eigen::Vector3d x1(m_cells_tet[0][0], m_cells_tet[0][1], m_cells_tet[0][2]);
            assert(is_step_valid(x0, x1));
        }
    } else if (m.top_simplex_type() == PrimitiveType::Triangle) {
        m_cells_tri.reserve(neighs.size());

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
            assert(vertices.size() == 3);

            std::array<double, 6> cell;
            for (size_t i = 0; i < 3; ++i) {
                const simplex::Simplex& v = vertices.simplex_vector()[i];
                const auto p = p_acc.const_vector_attribute(v.tuple());
                cell[2 * i + 0] = p[0];
                cell[2 * i + 1] = p[1];
            }
            m_cells_tri.emplace_back(cell);
        }

        // check if input is valid
        {
            Eigen::Vector2d x0(0, 0);
            Eigen::Vector2d x1(m_cells_tri[0][0], m_cells_tri[0][1]);
            assert(is_step_valid(x0, x1));
        }
    } else {
        log_and_throw_error("AmipsProblem unknown PrimitiveType");
    }
}

auto AmipsProblem::initial_value() const -> TVector
{
    if (m_pt == PrimitiveType::Tetrahedron) {
        Eigen::Vector3d tmp;
        for (int64_t d = 0; d < 3; ++d) {
            tmp(d) = m_cells_tet[0][d]; // m_cells[0] is the point we want to optimize
        }
        return tmp;
    } else if (m_pt == PrimitiveType::Triangle) {
        Eigen::Vector2d tmp;
        for (int64_t d = 0; d < 2; ++d) {
            tmp(d) = m_cells_tri[0][d]; // m_cells[0] is the point we want to optimize
        }
        return tmp;
    }
    log_and_throw_error("AmipsProblem unknown PrimitiveType");
}

double AmipsProblem::value(const TVector& x)
{
    assert(x.size() == 3 || x.size() == 2);

    double amips_energy = 0;
    if (m_pt == PrimitiveType::Tetrahedron) {
        for (auto c : m_cells_tet) {
            c[0] = x[0];
            c[1] = x[1];
            c[2] = x[2];
            amips_energy += wmtk::function::utils::Tet_AMIPS_energy(c);
        }
    } else if (m_pt == PrimitiveType::Triangle) {
        for (auto c : m_cells_tri) {
            c[0] = x[0];
            c[1] = x[1];
            amips_energy += wmtk::function::utils::Tri_AMIPS_energy(c);
        }
    }


    return amips_energy;
}

void AmipsProblem::gradient(const TVector& x, TVector& gradv)
{
    if (m_pt == PrimitiveType::Tetrahedron) {
        gradv.resize(3);
        gradv.setZero();
        Eigen::Vector3d tmp;

        Eigen::Vector3d amips_grad;
        amips_grad.setZero();
        for (auto c : m_cells_tet) {
            c[0] = x[0];
            c[1] = x[1];
            c[2] = x[2];
            wmtk::function::utils::Tet_AMIPS_jacobian(c, tmp);
            amips_grad += tmp;
        }

        gradv = amips_grad;
    } else if (m_pt == PrimitiveType::Triangle) {
        gradv.resize(2);
        gradv.setZero();
        Eigen::Vector2d tmp;

        Eigen::Vector2d amips_grad;
        amips_grad.setZero();
        for (auto c : m_cells_tri) {
            c[0] = x[0];
            c[1] = x[1];
            wmtk::function::utils::Tri_AMIPS_jacobian(c, tmp);
            amips_grad += tmp;
        }

        gradv = amips_grad;
    }
}

void AmipsProblem::hessian(const TVector& x, Eigen::MatrixXd& hessian)
{
    if (m_pt == PrimitiveType::Tetrahedron) {
        hessian.resize(3, 3);
        hessian.setZero();
        Eigen::Matrix3d tmp;

        Eigen::Matrix3d amips_hess;
        amips_hess.setZero();
        for (auto c : m_cells_tet) {
            c[0] = x[0];
            c[1] = x[1];
            c[2] = x[2];
            wmtk::function::utils::Tet_AMIPS_hessian(c, tmp);
            amips_hess += tmp;
        }

        hessian = amips_hess;
    } else if (m_pt == PrimitiveType::Triangle) {
        hessian.resize(2, 2);
        hessian.setZero();
        Eigen::Matrix2d tmp;

        Eigen::Matrix2d amips_hess;
        amips_hess.setZero();
        for (auto c : m_cells_tri) {
            c[0] = x[0];
            c[1] = x[1];
            wmtk::function::utils::Tri_AMIPS_hessian(c, tmp);
            amips_hess += tmp;
        }

        hessian = amips_hess;
    }
}

void AmipsProblem::solution_changed(const TVector& new_x) {}

bool AmipsProblem::is_step_valid(const TVector& x0, const TVector& x1)
{
    if (m_pt == PrimitiveType::Tetrahedron) {
        const Eigen::Vector3d p0 = x1;

        for (const auto& c : m_cells_tet) {
            const Eigen::Vector3d p1(c[3], c[4], c[5]);
            const Eigen::Vector3d p2(c[6], c[7], c[8]);
            const Eigen::Vector3d p3(c[9], c[10], c[11]);

            if (wmtk::utils::wmtk_orient3d(p3, p0, p1, p2) <= 0) {
                return false;
            }
        }

        return true;
    } else if (m_pt == PrimitiveType::Triangle) {
        const Eigen::Vector2d p0 = x1;

        for (const auto& c : m_cells_tri) {
            const Eigen::Vector2d p1(c[2], c[3]);
            const Eigen::Vector2d p2(c[4], c[5]);

            if (wmtk::utils::wmtk_orient2d(p0, p1, p2) <= 0) {
                return false;
            }
        }

        return true;
    }
    log_and_throw_error("AmipsProblem unknown PrimitiveType");
}
} // namespace wmtk::operations