#include "ConvergenceRate.hpp"

#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::utils {

ConvergenceRate::ConvergenceRate(
    const std::vector<double>& vals,
    const double relative_convergence_rate)
    : m_vals_curr(vals)
    , m_relative_convergence_rate(relative_convergence_rate)
{
    if (vals.empty()) {
        log_and_throw_error("Empty values in convergence rate are not allowed");
    }

    m_vals_prev = std::vector<double>(vals.size(), std::numeric_limits<double>::max());
}

void ConvergenceRate::update_vals(const std::vector<double>& vals)
{
    if (vals.size() != m_vals_curr.size()) {
        log_and_throw_error(
            "Cannot update values. Current values size: {}. New values size: {}",
            m_vals_curr.size(),
            vals.size());
    }

    m_vals_curr.swap(m_vals_prev);
    m_vals_curr = vals;
}

bool ConvergenceRate::is_converged()
{
    const auto rates = get_rates();

    for (const double r : rates) {
        if (r > m_relative_convergence_rate) {
            return false;
        }
    }

    return true;
}

void ConvergenceRate::print_rates()
{
    auto rates = get_rates();
    for (auto& r : rates) {
        r *= 100;
    }
    logger().info("Relative convergence rates: {:.4}%", fmt::join(rates, "%, "));
}

std::vector<double> ConvergenceRate::get_rates()
{
    std::vector<double> rates(m_vals_curr.size());

    for (size_t i = 0; i < m_vals_curr.size(); ++i) {
        rates[i] = m_vals_prev[i] - m_vals_curr[i];
    }

    return rates;
}

} // namespace wmtk::components::internal::utils
