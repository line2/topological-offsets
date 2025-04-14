#pragma once

#include <SimpleBVH/BVH.hpp>
#include <polysolve/nonlinear/Problem.hpp>
#include <wmtk/Mesh.hpp>

namespace wmtk::operations {

class AmipsOffsetProblem : public polysolve::nonlinear::Problem
{
public:
    using typename polysolve::nonlinear::Problem::Scalar;
    using typename polysolve::nonlinear::Problem::THessian;
    using typename polysolve::nonlinear::Problem::TVector;

    AmipsOffsetProblem(
        const simplex::Simplex& simplex,
        const attribute::MeshAttributeHandle& embedding_pos,
        const attribute::MeshAttributeHandle& offset_pos,
        const SimpleBVH::BVH& bvh,
        const double offset_distance);

    TVector initial_value() const;

    double value(const TVector& x) override;
    void gradient(const TVector& x, TVector& gradv) override;
    void hessian(const TVector& x, THessian& hessian) override
    {
        throw std::runtime_error("Sparse functions do not exist, use dense solver");
    }
    void hessian(const TVector& x, Eigen::MatrixXd& hessian) override;

    void solution_changed(const TVector& new_x) override;

    bool is_step_valid(const TVector& x0, const TVector& x1) override;

private:
    std::vector<std::array<double, 12>> m_cells;
    double m_inv_num_cells;

    Eigen::Vector3d m_p_opt;

    double m_w_amips = 1e-2;
    double m_w_offset = 1e4;
};

} // namespace wmtk::operations