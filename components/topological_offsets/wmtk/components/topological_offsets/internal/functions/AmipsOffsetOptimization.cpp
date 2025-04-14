#include "AmipsOffsetOptimization.hpp"

#include <Eigen/Core>
#include <polysolve/nonlinear/Solver.hpp>
#include <wmtk/Mesh.hpp>
#include <wmtk/Types.hpp>
#include <wmtk/attribute/Accessor.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/utils/Logger.hpp>

#include "AmipsOffsetProblem.hpp"


namespace wmtk::operations {

AmipsOffsetOptimization::AmipsOffsetOptimization(
    attribute::MeshAttributeHandle& embedding_pos,
    attribute::MeshAttributeHandle& offset_pos,
    const SimpleBVH::BVH& bvh,
    const double offset_distance)
    : AttributesUpdate(offset_pos.mesh())
    , m_embedding_pos(embedding_pos)
    , m_offset_pos(offset_pos)
    , mesh_embedding(embedding_pos.mesh())
    , mesh_offset(offset_pos.mesh())
    , m_bvh(bvh)
    , m_offset_distance(offset_distance)
    , m_amips(embedding_pos.mesh(), embedding_pos)
{
    assert(m_embedding_pos.holds<double>());
    assert(m_offset_pos.holds<double>());
    assert(mesh_embedding.top_simplex_type() == PrimitiveType::Tetrahedron);
    assert(mesh_offset.top_simplex_type() == PrimitiveType::Triangle);

    m_linear_solver_params = R"({"solver": "Eigen::LDLT"})"_json;
    m_nonlinear_solver_params = R"({"solver": "DenseNewton", "max_iterations": 10})"_json;

    create_solver();
}

void AmipsOffsetOptimization::create_solver()
{
    m_solver = polysolve::nonlinear::Solver::create(
        m_nonlinear_solver_params,
        m_linear_solver_params,
        1,
        opt_logger());
}


std::vector<simplex::Simplex> AmipsOffsetOptimization::execute(const simplex::Simplex& simplex)
{
    assert(mesh().top_simplex_type() == PrimitiveType::Triangle);

    AmipsOffsetProblem problem(simplex, m_embedding_pos, m_offset_pos, m_bvh, m_offset_distance);

    Eigen::VectorXd x = problem.initial_value();
    try {
        m_solver->minimize(problem, x);

        auto p_off_acc = mesh_offset.create_accessor(m_offset_pos.as<double>());
        auto p_emb_acc = mesh_embedding.create_accessor(m_embedding_pos.as<double>());

        p_off_acc.vector_attribute(simplex.tuple()) = x;
        const simplex::Simplex s_embedding = mesh_offset.map_to_parent(simplex);
        p_emb_acc.vector_attribute(s_embedding.tuple()) = x;

    } catch (const std::exception&) {
        return {};
    }

    return AttributesUpdate::execute(simplex);
}

} // namespace wmtk::operations