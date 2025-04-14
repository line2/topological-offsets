#include "AmipsOptimization.hpp"

#include <Eigen/Core>
#include <polysolve/nonlinear/Solver.hpp>
#include <wmtk/Mesh.hpp>
#include <wmtk/Types.hpp>
#include <wmtk/attribute/Accessor.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/utils/Logger.hpp>

#include "AmipsProblem.hpp"


namespace wmtk::operations {

AmipsOptimization::AmipsOptimization(attribute::MeshAttributeHandle& pos_handle)
    : AttributesUpdate(pos_handle.mesh())
    , m_pos_handle(pos_handle)
    , m_mesh(pos_handle.mesh())
    , m_amips(pos_handle.mesh(), pos_handle)
{
    assert(m_pos_handle.holds<double>());
    assert(
        m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron ||
        m_mesh.top_simplex_type() == PrimitiveType::Triangle);

    m_linear_solver_params = R"({"solver": "Eigen::LDLT"})"_json;
    m_nonlinear_solver_params = R"({"solver": "DenseNewton", "max_iterations": 10})"_json;

    create_solver();
}

void AmipsOptimization::create_solver()
{
    m_solver = polysolve::nonlinear::Solver::create(
        m_nonlinear_solver_params,
        m_linear_solver_params,
        1,
        opt_logger());
}


std::vector<simplex::Simplex> AmipsOptimization::execute(const simplex::Simplex& simplex)
{
    assert(
        mesh().top_simplex_type() == PrimitiveType::Tetrahedron ||
        m_mesh.top_simplex_type() == PrimitiveType::Triangle);

    AmipsProblem problem(simplex, m_pos_handle);

    Eigen::VectorXd x = problem.initial_value();
    try {
        m_solver->minimize(problem, x);

    } catch (const std::exception&) {
        // PolySolve might fail but still returns a valid position
        // return {};
    }

    auto p_acc = m_mesh.create_accessor(m_pos_handle.as<double>());
    p_acc.vector_attribute(simplex.tuple()) = x;

    return AttributesUpdate::execute(simplex);
}

} // namespace wmtk::operations