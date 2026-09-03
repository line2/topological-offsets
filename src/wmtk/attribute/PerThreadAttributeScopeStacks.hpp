#pragma once

#include <tbb/combinable.h>

#include "internal/AttributeTransactionStack.hpp"
namespace wmtk::attribute {
template <typename T>
class PerThreadAttributeScopeStacks
{
public:
    PerThreadAttributeScopeStacks() = default;
    PerThreadAttributeScopeStacks(PerThreadAttributeScopeStacks&&) = default;
    PerThreadAttributeScopeStacks& operator=(PerThreadAttributeScopeStacks&&) = default;
    internal::AttributeTransactionStack<T>& local();
    const internal::AttributeTransactionStack<T>& local() const;

private:
    // one transaction stack per thread so that concurrent operations keep
    // their scope bookkeeping isolated
    mutable tbb::combinable<internal::AttributeTransactionStack<T>> m_stacks;
};


template <typename T>
inline internal::AttributeTransactionStack<T>& PerThreadAttributeScopeStacks<T>::local()
{
    return m_stacks.local();
}
template <typename T>
inline const internal::AttributeTransactionStack<T>& PerThreadAttributeScopeStacks<T>::local() const
{
    return m_stacks.local();
}
} // namespace wmtk::attribute
