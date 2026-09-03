#pragma once

#include <wmtk/utils/Rational.hpp>
#include "AttributeScope.hpp"
#include "internal/AttributeTransactionStack.hpp"
#include "CachingAccessor.hpp"

namespace wmtk::attribute {

template <typename T, int Dim>
inline CachingAccessor<T, Dim>::CachingAccessor(
    Mesh& mesh_in,
    const TypedAttributeHandle<T>& handle)
    : BaseType(mesh_in, handle)
{}
template <typename T, int Dim>
CachingAccessor<T, Dim>::CachingAccessor(const Mesh& mesh_in, const TypedAttributeHandle<T>& handle)
    : BaseType(mesh_in, handle)
{}

template <typename T, int Dim>
inline CachingAccessor<T, Dim>::~CachingAccessor() = default;

template <typename T, int Dim>
inline bool CachingAccessor<T, Dim>::has_stack() const
{
    return !cache_stack().empty();
}

template <typename T, int Dim>
inline bool CachingAccessor<T, Dim>::writing_enabled() const
{
    return cache_stack().writing_enabled();
}

template <typename T, int Dim>
inline int64_t CachingAccessor<T, Dim>::stack_depth() const
{
    return cache_stack().size();
}

template <typename T, int Dim>
template <int D>
inline auto CachingAccessor<T, Dim>::vector_attribute(const int64_t index) -> MapResult<D>
{
    return cache_stack().template vector_attribute<D>(*this, index);
    // return BaseType::template vector_attribute<D>( index);
}


template <typename T, int Dim>
inline auto CachingAccessor<T, Dim>::scalar_attribute(const int64_t index) -> T&
{
    return cache_stack().scalar_attribute(*this, index);
}

template <typename T, int Dim>
template <int D>
inline auto CachingAccessor<T, Dim>::const_vector_attribute(const int64_t index) const
    -> ConstMapResult<D>
{
    return cache_stack().template const_vector_attribute<D>(*this, index);
}


template <typename T, int Dim>
inline auto CachingAccessor<T, Dim>::const_scalar_attribute(const int64_t index) const -> T
{
    return cache_stack().const_scalar_attribute(*this, index);
}

template <typename T, int Dim>
inline auto CachingAccessor<T, Dim>::vector_attribute(const int64_t index) const -> ConstMapResult<>
{
    return const_vector_attribute(index);
}
template <typename T, int Dim>
inline T CachingAccessor<T, Dim>::scalar_attribute(const int64_t index) const
{
    return const_scalar_attribute(index);
}

template <typename T, int Dim>
inline auto CachingAccessor<T, Dim>::scalar_attribute(const int64_t index, int8_t offset) -> T&
{
    return cache_stack().scalar_attribute(*this, index, offset);
}


template <typename T, int Dim>
inline auto CachingAccessor<T, Dim>::const_scalar_attribute(const int64_t index, int8_t offset)
    const -> T
{
    return cache_stack().const_scalar_attribute(*this, index, offset);
}

// template class CachingAccessor<char>;
// template class CachingAccessor<int64_t>;
// template class CachingAccessor<double>;
// template class CachingAccessor<Rational>;
} // namespace wmtk::attribute
