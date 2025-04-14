#pragma once

#include <wmtk/attribute/TypedAttributeHandle.hpp>
#include "Invariant.hpp"

namespace wmtk {
template <typename AttributeType>
class Swap23EnergyBeforeInvariant : public Invariant
{
public:
    Swap23EnergyBeforeInvariant(
        const Mesh& m,
        const attribute::TypedAttributeHandle<AttributeType>& coordinate);
    using Invariant::Invariant;

    bool before(const simplex::Simplex& t) const override;

private:
    const attribute::TypedAttributeHandle<AttributeType> m_coordinate_handle;
};
} // namespace wmtk
