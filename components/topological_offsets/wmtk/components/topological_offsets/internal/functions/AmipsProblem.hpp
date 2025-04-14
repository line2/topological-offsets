#pragma once

#include <SimpleBVH/BVH.hpp>
#include <polysolve/nonlinear/Problem.hpp>
#include <wmtk/Mesh.hpp>

namespace wmtk::operations {

class AmipsProblem : public polysolve::nonlinear::Problem
{
public:
    using typename polysolve::nonlinear::Problem::Scalar;
    using typename polysolve::nonlinear::Problem::THessian;
    using typename polysolve::nonlinear::Problem::TVector;

    AmipsProblem(const simplex::Simplex& simplex, const attribute::MeshAttributeHandle& pos_handle);

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
    std::vector<std::array<double, 12>> m_cells_tet;
    std::vector<std::array<double, 6>> m_cells_tri;
    PrimitiveType m_pt;
};

} // namespace wmtk::operations