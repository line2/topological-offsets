#pragma once

#include <SimpleBVH/BVH.hpp>
#include <polysolve/Types.hpp>
#include <tbb/combinable.h>
#include <wmtk/function/Function.hpp>
#include <wmtk/function/simplex/AMIPS.hpp>
#include <wmtk/operations/AttributesUpdate.hpp>


namespace polysolve::nonlinear {
class Solver;
}

namespace wmtk::function {
class Function;
}

namespace wmtk::operations {

class AmipsOptimization : public AttributesUpdate
{
public:
    AmipsOptimization(attribute::MeshAttributeHandle& pos_handle);

    std::vector<simplex::Simplex> execute(const simplex::Simplex& simplex) override;

    const polysolve::json& linear_solver_params() const { return m_linear_solver_params; }
    const polysolve::json& nonlinear_solver_params() const { return m_nonlinear_solver_params; }


    void set_linear_solver_params(const polysolve::json& params)
    {
        m_linear_solver_params = params;
        create_solver();
    }

    void set_nonlinear_solver_params(const polysolve::json& params)
    {
        m_nonlinear_solver_params = params;
        create_solver();
    }

private:
    // One solver per thread (solvers carry internal state and are not
    // thread-safe); created lazily on first use of each thread.
    tbb::combinable<std::shared_ptr<polysolve::nonlinear::Solver>> m_thread_solvers;
    std::shared_ptr<polysolve::nonlinear::Solver> thread_solver();
    const attribute::MeshAttributeHandle& m_pos_handle;
    Mesh& m_mesh;

    polysolve::json m_linear_solver_params;
    polysolve::json m_nonlinear_solver_params;

    function::AMIPS m_amips;

    void create_solver();
};

} // namespace wmtk::operations