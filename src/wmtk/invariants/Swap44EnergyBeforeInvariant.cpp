#include "Swap44EnergyBeforeInvariant.hpp"
#include <wmtk/Mesh.hpp>
#include <wmtk/function/utils/amips.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/utils/orient.hpp>

namespace wmtk {
template <typename AttributeType>
Swap44EnergyBeforeInvariant<AttributeType>::Swap44EnergyBeforeInvariant(
    const Mesh& m,
    const attribute::TypedAttributeHandle<AttributeType>& coordinate,
    int64_t collapse_index,
    double eps)
    : Invariant(m, true, false, false)
    , m_coordinate_handle(coordinate)
    , m_collapse_index(collapse_index)
    , m_eps(eps)
{}

template <typename AttributeType>
bool Swap44EnergyBeforeInvariant<AttributeType>::before(const simplex::Simplex& t) const
{
    assert(simplex::top_dimension_cofaces(mesh(), t).size() == 4);
    constexpr static PrimitiveType PV = PrimitiveType::Vertex;
    constexpr static PrimitiveType PE = PrimitiveType::Edge;
    constexpr static PrimitiveType PF = PrimitiveType::Triangle;
    constexpr static PrimitiveType PT = PrimitiveType::Tetrahedron;

    auto accessor = mesh().create_const_accessor(m_coordinate_handle);

    // get the coords of the vertices
    // input edge end points
    const Tuple e0 = t.tuple();
    const Tuple e1 = mesh().switch_tuple(e0, PV);
    // other four vertices
    std::array<Tuple, 4> v;
    auto iter_tuple = e0;
    for (int64_t i = 0; i < 4; ++i) {
        v[i] = mesh().switch_tuples(iter_tuple, {PE, PV});
        iter_tuple = mesh().switch_tuples(iter_tuple, {PF, PT});
    }
    assert(iter_tuple == e0);

    std::array<Eigen::Vector3<AttributeType>, 6> positions = {
        {accessor.const_vector_attribute(v[(m_collapse_index + 0) % 4]),
         accessor.const_vector_attribute(v[(m_collapse_index + 1) % 4]),
         accessor.const_vector_attribute(v[(m_collapse_index + 2) % 4]),
         accessor.const_vector_attribute(v[(m_collapse_index + 3) % 4]),
         accessor.const_vector_attribute(e0),
         accessor.const_vector_attribute(e1)}};

    std::array<Eigen::Vector3d, 6> positions_double;
    if constexpr (std::is_same<AttributeType, Rational>::value) {
        positions_double = {
            {positions[0].template cast<double>(),
             positions[1].template cast<double>(),
             positions[2].template cast<double>(),
             positions[3].template cast<double>(),
             positions[4].template cast<double>(),
             positions[5].template cast<double>()}};
    } else {
        positions_double = positions;
    }

    std::array<std::array<int, 4>, 4> old_tets = {
        {{{0, 1, 4, 5}}, {{1, 2, 4, 5}}, {{2, 3, 4, 5}}, {{3, 0, 4, 5}}}};
    std::array<std::array<int, 4>, 4> new_tets = {
        {{{0, 1, 2, 4}}, {{0, 2, 3, 4}}, {{0, 1, 2, 5}}, {{0, 2, 3, 5}}}};

    double old_energy_max = std::numeric_limits<double>::lowest();
    double new_energy_max = std::numeric_limits<double>::lowest();

    for (int i = 0; i < 4; ++i) {
        if (utils::wmtk_orient3d(
                positions[old_tets[i][0]],
                positions[old_tets[i][1]],
                positions[old_tets[i][2]],
                positions[old_tets[i][3]]) > 0) {
            auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                positions_double[old_tets[i][0]][0],
                positions_double[old_tets[i][0]][1],
                positions_double[old_tets[i][0]][2],
                positions_double[old_tets[i][1]][0],
                positions_double[old_tets[i][1]][1],
                positions_double[old_tets[i][1]][2],
                positions_double[old_tets[i][2]][0],
                positions_double[old_tets[i][2]][1],
                positions_double[old_tets[i][2]][2],
                positions_double[old_tets[i][3]][0],
                positions_double[old_tets[i][3]][1],
                positions_double[old_tets[i][3]][2],
            }});

            old_energy_max = std::max(energy, old_energy_max);
        } else {
            auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                positions_double[old_tets[i][1]][0],
                positions_double[old_tets[i][1]][1],
                positions_double[old_tets[i][1]][2],
                positions_double[old_tets[i][0]][0],
                positions_double[old_tets[i][0]][1],
                positions_double[old_tets[i][0]][2],
                positions_double[old_tets[i][2]][0],
                positions_double[old_tets[i][2]][1],
                positions_double[old_tets[i][2]][2],
                positions_double[old_tets[i][3]][0],
                positions_double[old_tets[i][3]][1],
                positions_double[old_tets[i][3]][2],
            }});

            old_energy_max = std::max(energy, old_energy_max);
        }

        if (utils::wmtk_orient3d(
                positions[new_tets[i][0]],
                positions[new_tets[i][1]],
                positions[new_tets[i][2]],
                positions[new_tets[i][3]]) > 0) {
            auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                positions_double[new_tets[i][0]][0],
                positions_double[new_tets[i][0]][1],
                positions_double[new_tets[i][0]][2],
                positions_double[new_tets[i][1]][0],
                positions_double[new_tets[i][1]][1],
                positions_double[new_tets[i][1]][2],
                positions_double[new_tets[i][2]][0],
                positions_double[new_tets[i][2]][1],
                positions_double[new_tets[i][2]][2],
                positions_double[new_tets[i][3]][0],
                positions_double[new_tets[i][3]][1],
                positions_double[new_tets[i][3]][2],
            }});

            new_energy_max = std::max(energy, new_energy_max);
        } else {
            auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                positions_double[new_tets[i][1]][0],
                positions_double[new_tets[i][1]][1],
                positions_double[new_tets[i][1]][2],
                positions_double[new_tets[i][0]][0],
                positions_double[new_tets[i][0]][1],
                positions_double[new_tets[i][0]][2],
                positions_double[new_tets[i][2]][0],
                positions_double[new_tets[i][2]][1],
                positions_double[new_tets[i][2]][2],
                positions_double[new_tets[i][3]][0],
                positions_double[new_tets[i][3]][1],
                positions_double[new_tets[i][3]][2],
            }});

            new_energy_max = std::max(energy, new_energy_max);
        }
    }

    return old_energy_max > new_energy_max * m_eps;
}

template class Swap44EnergyBeforeInvariant<Rational>;
template class Swap44EnergyBeforeInvariant<double>;

} // namespace wmtk
