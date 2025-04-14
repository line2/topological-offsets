#pragma once

#include <vector>

namespace wmtk::components::internal::utils {
class ConvergenceRate
{
public:
    ConvergenceRate(const std::vector<double>& vals, const double relative_convergence_rate);

    void update_vals(const std::vector<double>& vals);

    bool is_converged();

    void print_rates();

private:
    std::vector<double> get_rates();

private:
    std::vector<double> m_vals_curr; // current
    std::vector<double> m_vals_prev; // previous
    double m_relative_convergence_rate;
};


} // namespace wmtk::components::internal::utils
