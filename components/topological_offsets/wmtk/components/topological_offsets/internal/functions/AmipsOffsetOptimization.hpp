#pragma once

#include <SimpleBVH/BVH.hpp>
#include <polysolve/Types.hpp>
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

class AmipsOffsetOptimization : public AttributesUpdate
{
public:
    AmipsOffsetOptimization(
        attribute::MeshAttributeHandle& embedding_pos,
        attribute::MeshAttributeHandle& offset_pos,
        const SimpleBVH::BVH& bvh,
        const double offset_distance);

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
    std::shared_ptr<polysolve::nonlinear::Solver> m_solver;
    const attribute::MeshAttributeHandle& m_embedding_pos;
    const attribute::MeshAttributeHandle& m_offset_pos;
    Mesh& mesh_embedding;
    Mesh& mesh_offset;
    const SimpleBVH::BVH& m_bvh;
    const double m_offset_distance;

    polysolve::json m_linear_solver_params;
    polysolve::json m_nonlinear_solver_params;

    function::AMIPS m_amips;

    void create_solver();
};

} // namespace wmtk::operations