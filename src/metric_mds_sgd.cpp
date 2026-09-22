// Pair updates adapted from Zheng et al. (2018) and Hangan et al. (2026).
// Source revisions and deliberate differences: dev/design/metric-mds-sgd.md.
// Third-party MIT and BSD-3-Clause notices: inst/COPYRIGHTS.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

namespace {
struct Pair { int i, j; double target; };
struct Score { double raw, profiled; };

std::size_t checked_product(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max() / b)
        Rcpp::stop("SGD workspace size overflow");
    return a * b;
}

// Standardized engine plus explicit bounded sampling avoids implementation-
// dependent std::shuffle/uniform_int_distribution behavior.
std::uint64_t bounded(std::mt19937_64& rng, std::uint64_t bound) {
    const std::uint64_t threshold = -bound % bound;
    std::uint64_t x;
    do { x = rng(); } while (x < threshold);
    return x % bound;
}

double separation(const Rcpp::NumericMatrix& x, int i, int j) {
    double distance = 0;
    for (int d = 0; d < x.ncol(); ++d)
        distance = std::hypot(distance, x(i, d) - x(j, d));
    if (!std::isfinite(distance)) Rcpp::stop("Nonfinite SGD separation");
    return distance;
}

Score score(const Rcpp::NumericMatrix& x, const std::vector<Pair>& pairs) {
    long double raw = 0, cross = 0, squares = 0;
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
        const Pair& p = pairs[k];
        const long double r = separation(x, p.i, p.j);
        const long double error = r - p.target;
        raw += error * error;
        cross += r * p.target;
        squares += r * r;
    }
    if (!(squares > 0)) Rcpp::stop("Collapsed SGD checkpoint");
    const long double scale = cross / squares;
    long double profiled = 0;
    // Sum residuals directly rather than subtracting nearly equal large sums.
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
        const Pair& p = pairs[k];
        const long double error = scale * separation(x, p.i, p.j) - p.target;
        profiled += error * error;
    }
    const Score result{static_cast<double>(raw), static_cast<double>(profiled)};
    if (!std::isfinite(result.raw) || !std::isfinite(result.profiled))
        Rcpp::stop("Nonfinite SGD checkpoint stress");
    return result;
}
}

// Internal Rcpp entry point; shuffle=false exists for deterministic kernel tests.
// [[Rcpp::export]]
Rcpp::List grip_sgd_mds_cpp(Rcpp::NumericMatrix start,
                          Rcpp::NumericVector targets,
                          Rcpp::NumericVector rates, int seed,
                          int checkpoint_every, double max_workspace_bytes,
                          bool shuffle = true) {
    const int n = start.nrow(), dim = start.ncol();
    if (n < 2 || dim < 1 || checkpoint_every < 1 || rates.size() < 1 || seed < 0 ||
        !std::isfinite(max_workspace_bytes) || max_workspace_bytes <= 0)
        Rcpp::stop("Invalid native SGD controls or dimensions");
    if (rates.size() > std::numeric_limits<int>::max()) Rcpp::stop("Too many SGD epochs");
    const std::size_t count = checked_product(static_cast<std::size_t>(n), n - 1) / 2;
    if (static_cast<std::uint64_t>(targets.size()) != count)
        Rcpp::stop("SGD target count does not match the configuration");
    const int epochs = static_cast<int>(rates.size());
    const std::size_t checkpoints = 1 + static_cast<std::size_t>(epochs / checkpoint_every) +
        (epochs % checkpoint_every != 0);
    const long double workspace = static_cast<long double>(checked_product(count, sizeof(Pair) + sizeof(std::size_t))) +
        2.0L * checked_product(static_cast<std::size_t>(start.size()), sizeof(double)) +
        checked_product(checkpoints, 5 * sizeof(double)) + dim * sizeof(double);
    if (workspace > max_workspace_bytes)
        Rcpp::stop("SGD native workspace exceeds max_workspace_bytes; R input matrices require additional memory");
    for (R_xlen_t k = 0; k < start.size(); ++k) {
        if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
        if (!std::isfinite(start[k])) Rcpp::stop("Nonfinite SGD start");
    }
    for (double rate : rates)
        if (!std::isfinite(rate) || rate <= 0) Rcpp::stop("SGD rates must be finite and positive");
    std::vector<Pair> pairs;
    pairs.reserve(count);
    R_xlen_t k = 0;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j, ++k) {
            if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
            const double target = targets[k];
            if (!std::isfinite(target) || target < 0) Rcpp::stop("Invalid SGD target distance");
            pairs.push_back(Pair{i, j, target});
        }
    }
    Rcpp::NumericMatrix x = Rcpp::clone(start), best = Rcpp::clone(start);
    Score initial = score(x, pairs);
    double best_loss = initial.profiled;
    int best_epoch = 0;
    std::vector<std::size_t> order(count);
    std::vector<double> direction(dim);
    std::mt19937_64 rng(static_cast<std::uint64_t>(seed));
    Rcpp::NumericVector epoch_trace(checkpoints), raw_trace(checkpoints),
        profiled_trace(checkpoints), updates_trace(checkpoints), rate_trace(checkpoints);
    raw_trace[0] = initial.raw;
    profiled_trace[0] = initial.profiled;
    rate_trace[0] = NA_REAL;
    std::size_t row = 1;
    double updates = 0;
    for (int epoch = 0; epoch < epochs; ++epoch) {
        Rcpp::checkUserInterrupt();
        std::iota(order.begin(), order.end(), std::size_t(0));
        if (shuffle) {
            for (std::size_t remaining = count; remaining > 1; --remaining) {
                if ((remaining & 16383) == 0) Rcpp::checkUserInterrupt();
                std::swap(order[remaining - 1], order[bounded(rng, remaining)]);
            }
        }
        const double mu = std::min(rates[epoch], 1.0);
        for (std::size_t pos = 0; pos < count; ++pos) {
            if ((pos & 16383) == 0) Rcpp::checkUserInterrupt();
            const Pair& p = pairs[order[pos]];
            const double distance = separation(x, p.i, p.j);
            if (distance == 0 && p.target == 0) continue;
            if (distance == 0) {
                double norm = 0;
                do {
                    norm = 0;
                    for (int d = 0; d < dim; ++d) {
                        direction[d] = 2.0 * static_cast<double>(rng() >> 11) * 0x1.0p-53 - 1.0;
                        norm = std::hypot(norm, direction[d]);
                    }
                } while (norm == 0);
                for (int d = 0; d < dim; ++d) direction[d] /= norm;
            } else {
                for (int d = 0; d < dim; ++d)
                    direction[d] = (x(p.i, d) - x(p.j, d)) / distance;
            }
            const double move = (distance - p.target) * (0.5 * mu);
            for (int d = 0; d < dim; ++d) {
                const double change = move * direction[d];
                x(p.i, d) -= change;
                x(p.j, d) += change;
                if (!std::isfinite(x(p.i, d)) || !std::isfinite(x(p.j, d)))
                    Rcpp::stop("Nonfinite SGD coordinate update");
            }
        }
        updates += static_cast<double>(count); // Includes examined zero/zero pairs.
        const int completed = epoch + 1;
        if (completed % checkpoint_every == 0 || completed == epochs) {
            const Score current = score(x, pairs);
            epoch_trace[row] = completed;
            raw_trace[row] = current.raw;
            profiled_trace[row] = current.profiled;
            updates_trace[row] = updates;
            rate_trace[row++] = rates[epoch];
            if (current.profiled < best_loss) {
                std::copy(x.begin(), x.end(), best.begin());
                best_loss = current.profiled;
                best_epoch = completed;
            }
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("conf") = best,
        Rcpp::Named("terminal_conf") = x,
        Rcpp::Named("niter") = epochs,
        Rcpp::Named("stress") = best_loss,
        Rcpp::Named("best_epoch") = best_epoch,
        Rcpp::Named("pair_updates") = updates,
        Rcpp::Named("native_seed") = seed,
        Rcpp::Named("workspace_bytes") = static_cast<double>(workspace),
        Rcpp::Named("trace") = Rcpp::DataFrame::create(
            Rcpp::Named("epoch") = epoch_trace,
            Rcpp::Named("raw_stress") = raw_trace,
            Rcpp::Named("profiled_stress") = profiled_trace,
            Rcpp::Named("pair_updates") = updates_trace,
            Rcpp::Named("learning_rate") = rate_trace));
}
