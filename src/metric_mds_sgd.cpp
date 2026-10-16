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
#include <sstream>
#include <iomanip>

namespace {
struct Pair { int i, j; double target, weight, reverse; };
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

Score score(const Rcpp::NumericMatrix& x, const std::vector<Pair>& pairs, bool profile = true) {
    long double raw = 0, cross = 0, squares = 0;
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
        const Pair& p = pairs[k];
        const long double r = separation(x, p.i, p.j);
        const long double error = r - p.target;
        const double weight = p.weight / 2 + p.reverse / 2;
        raw += weight * error * error;
        cross += weight * r * p.target;
        squares += weight * r * r;
    }
    if (!profile) {
        if (!std::isfinite(static_cast<double>(raw))) Rcpp::stop("Nonfinite SGD checkpoint stress");
        return Score{static_cast<double>(raw), NA_REAL};
    }
    if (!(squares > 0)) Rcpp::stop("Collapsed SGD checkpoint");
    const long double scale = cross / squares;
    long double profiled = 0;
    // Sum residuals directly rather than subtracting nearly equal large sums.
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
        const Pair& p = pairs[k];
        const long double error = scale * separation(x, p.i, p.j) - p.target;
        profiled += (p.weight / 2 + p.reverse / 2) * error * error;
    }
    const Score result{static_cast<double>(raw), static_cast<double>(profiled)};
    if (!std::isfinite(result.raw) || !std::isfinite(result.profiled))
        Rcpp::stop("Nonfinite SGD checkpoint stress");
    return result;
}
// Weighted mean update field, independent of the current learning rate.
// This is a practical stationarity diagnostic, not a gradient certificate for
// asymmetric endpoint weights. A positive-target collision cannot pass it.
double update_residual(const Rcpp::NumericMatrix& x, const std::vector<Pair>& pairs) {
    const int n=x.nrow(), dim=x.ncol();
    std::vector<long double> total(n,0), force(static_cast<std::size_t>(n)*dim,0);
    for (const auto& p : pairs) { total[p.i]+=p.weight; total[p.j]+=p.reverse; }
    for (std::size_t k=0;k<pairs.size();++k) {
        if ((k & 16383)==0) Rcpp::checkUserInterrupt();
        const auto& p=pairs[k]; const double d=separation(x,p.i,p.j);
        if (d==0) { if (p.target>0) return R_PosInf; else continue; }
        for (int j=0;j<dim;++j) {
            const long double r=(d-p.target)*((x(p.i,j)-x(p.j,j))/d);
            if (total[p.i]>0) force[p.i*dim+j]-=(p.weight/total[p.i])*r;
            if (total[p.j]>0) force[p.j*dim+j]+=(p.reverse/total[p.j])*r;
        }
    }
    double maximum=0;
    for (int i=0;i<n;++i) {
        double norm=0;
        for (int j=0;j<dim;++j) norm=std::hypot(norm,static_cast<double>(force[i*dim+j]));
        maximum=std::max(maximum,norm);
    }
    return maximum;
}

double centered_motion(const Rcpp::NumericMatrix& x,const Rcpp::NumericMatrix& previous) {
    long double sum=0;
    for(int j=0;j<x.ncol();++j) {
        long double mean=0;
        for(int i=0;i<x.nrow();++i) mean+=x(i,j)-previous(i,j);
        mean/=x.nrow();
        for(int i=0;i<x.nrow();++i) {
            const long double d=x(i,j)-previous(i,j)-mean; sum+=d*d;
        }
    }
    return std::sqrt(static_cast<double>(sum/x.nrow()));
}

}

// Internal Rcpp entry point; shuffle=false exists for deterministic kernel tests.
// [[Rcpp::export]]
Rcpp::List grip_sgd_mds_cpp(Rcpp::NumericMatrix start,
                          Rcpp::NumericVector targets,
                          Rcpp::NumericVector rates, int seed,
                          int checkpoint_every, double max_workspace_bytes,
                          bool shuffle = true,
                          Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue,
                          Rcpp::Nullable<Rcpp::IntegerMatrix> endpoints = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericVector> reverse_weights = R_NilValue,
                          bool retain_best = true,
                          Rcpp::Nullable<Rcpp::List> state = R_NilValue,
                          Rcpp::Nullable<Rcpp::List> stopping = R_NilValue,
                          int schedule_epochs = 0, bool keep_state = false) {
    const int n = start.nrow(), dim = start.ncol();
    if (n < 2 || dim < 1 || checkpoint_every < 1 || rates.size() < 1 || seed < 0 ||
        !std::isfinite(max_workspace_bytes) || max_workspace_bytes <= 0)
        Rcpp::stop("Invalid native SGD controls or dimensions");
    if (rates.size() > std::numeric_limits<int>::max()) Rcpp::stop("Too many SGD epochs");
    const bool sparse = endpoints.isNotNull();
    const Rcpp::IntegerMatrix indices = sparse ? Rcpp::IntegerMatrix(endpoints) : Rcpp::IntegerMatrix(0,2);
    if (indices.ncol() != 2) Rcpp::stop("SGD endpoints must have two columns");
    const std::size_t count = sparse ? indices.nrow() : checked_product(static_cast<std::size_t>(n), n - 1) / 2;
    if (!count) Rcpp::stop("SGD needs at least one pair");
    if (static_cast<std::uint64_t>(targets.size()) != count)
        Rcpp::stop("SGD target count does not match the configuration");
    const bool weighted = weights.isNotNull();
    const Rcpp::NumericVector pair_weights = weighted ? Rcpp::NumericVector(weights) : Rcpp::NumericVector();
    if (weighted && static_cast<std::uint64_t>(pair_weights.size()) != count)
        Rcpp::stop("SGD weight count does not match the configuration");
    const Rcpp::NumericVector reverse = reverse_weights.isNotNull() ? Rcpp::NumericVector(reverse_weights) : pair_weights;
    if (reverse_weights.isNotNull() && (!weighted || static_cast<std::uint64_t>(reverse.size()) != count))
        Rcpp::stop("SGD reverse weight count does not match the configuration");
    const int epochs = static_cast<int>(rates.size());
    Rcpp::List saved = state.isNotNull() ? Rcpp::List(state) : Rcpp::List();
    const int offset = state.isNotNull() ? Rcpp::as<int>(saved["epoch"]) : 0;
    if (offset<0 || offset>std::numeric_limits<int>::max()-epochs)
        Rcpp::stop("Invalid continuation epoch range");
    if (schedule_epochs<=0) schedule_epochs=offset+epochs;
    const bool assess=stopping.isNotNull();
    Rcpp::List rule=assess ? Rcpp::List(stopping) : Rcpp::List();
    const double rel_tol=assess ? Rcpp::as<double>(rule["relative.tol"]) : 0;
    const double motion_tol=assess ? Rcpp::as<double>(rule["motion.tol"]) : 0;
    const double residual_tol=assess ? Rcpp::as<double>(rule["residual.tol"]) : 0;
    const int patience=assess ? Rcpp::as<int>(rule["patience"]) : 1;
    const int min_epochs=assess ? Rcpp::as<int>(rule["min.epochs"]) : 1;
    if (assess && (!std::isfinite(rel_tol)||rel_tol<0 || !std::isfinite(motion_tol)||motion_tol<0 ||
        !std::isfinite(residual_tol)||residual_tol<0 || patience<1 || min_epochs<1))
        Rcpp::stop("Invalid native stopping rule");
    // Two extra rows cover a schedule endpoint and a non-checkpoint call boundary.
    const std::size_t checkpoints = 3 + static_cast<std::size_t>(epochs / checkpoint_every);
    const long double workspace = static_cast<long double>(checked_product(count, sizeof(Pair) + sizeof(std::size_t))) +
        (assess ? 5.0L : 3.0L) * checked_product(static_cast<std::size_t>(start.size()), sizeof(double)) +
        checked_product(checkpoints, 10 * sizeof(double)) + dim * sizeof(double) +
        (assess ? static_cast<long double>(n)*(dim+1)*sizeof(long double) : 0);
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
    int dense_i = 0, dense_j = 1;
    for (std::size_t k=0; k<count; ++k) {
        if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
        if (sparse && (indices(k,0) == NA_INTEGER || indices(k,1) == NA_INTEGER))
            Rcpp::stop("Invalid SGD endpoint pair");
        const int i = sparse ? indices(k,0)-1 : dense_i;
        const int j = sparse ? indices(k,1)-1 : dense_j;
        if (i < 0 || j <= i || j >= n) Rcpp::stop("Invalid SGD endpoint pair");
        const double target = targets[k];
        if (!std::isfinite(target) || target < 0) Rcpp::stop("Invalid SGD target distance");
        const double weight = weighted ? pair_weights[k] : 1.0;
        const double back = weighted ? reverse[k] : 1.0;
        if (!std::isfinite(weight) || !std::isfinite(back) || weight < 0 || back < 0 ||
            (!sparse && (weight == 0 || back == 0)) || (weight == 0 && back == 0))
            Rcpp::stop("SGD pair weights must be finite and positive at an active endpoint");
        pairs.push_back(Pair{i,j,target,weight,back});
        if (!sparse && ++dense_j == n) { ++dense_i; dense_j = dense_i+1; }
    }
    Rcpp::NumericMatrix x = Rcpp::clone(start), best = Rcpp::clone(start);
    Score initial = score(x, pairs, retain_best);
    double best_loss = retain_best ? initial.profiled : initial.raw;
    int best_epoch = 0, stable_count=0, assessments=0;
    double previous_score=best_loss;
    Rcpp::NumericMatrix previous = assess ? Rcpp::clone(x) : Rcpp::NumericMatrix(0,0);
    std::vector<std::size_t> order(count);
    std::vector<double> direction(dim);
    std::mt19937_64 rng(static_cast<std::uint64_t>(seed));
    if (state.isNotNull()) {
        Rcpp::NumericMatrix old_best=saved["best"];
        if(old_best.nrow()!=n || old_best.ncol()!=dim) Rcpp::stop("Invalid continuation shape");
        best=Rcpp::clone(old_best); best_loss=Rcpp::as<double>(saved["best_loss"]);
        best_epoch=Rcpp::as<int>(saved["best_epoch"]);
        stable_count=Rcpp::as<int>(saved["stable_count"]);
        assessments=Rcpp::as<int>(saved["assessments"]);
        previous_score=Rcpp::as<double>(saved["checkpoint_score"]);
        if(assess) {
            Rcpp::NumericMatrix old_previous=saved["checkpoint_coords"];
            if(old_previous.nrow()!=n || old_previous.ncol()!=dim) Rcpp::stop("Invalid stopping state shape");
            previous=Rcpp::clone(old_previous);
        }
        std::istringstream stream(Rcpp::as<std::string>(saved["rng"]));
        if (!(stream>>rng)) Rcpp::stop("Invalid continuation RNG state");
    }
    std::vector<double> epoch_trace,raw_trace,profiled_trace,updates_trace,rate_trace,
        motion_trace,residual_trace,relative_trace;
    std::vector<int> scheduled_trace,streak_trace;
    auto record = [&](int epoch, Score value, double rate, bool scheduled,
                      double motion, double residual, double relative) {
        epoch_trace.push_back(epoch); raw_trace.push_back(value.raw); profiled_trace.push_back(value.profiled);
        updates_trace.push_back(static_cast<double>(epoch)*count);rate_trace.push_back(rate);
        scheduled_trace.push_back(scheduled);motion_trace.push_back(motion);residual_trace.push_back(residual);
        relative_trace.push_back(relative);streak_trace.push_back(stable_count);
    };
    if (!state.isNotNull()) record(0,initial,NA_REAL,true,NA_REAL,NA_REAL,NA_REAL);
    double updates=static_cast<double>(offset)*count;
    int completed=offset;bool stabilized=false;
    for (int epoch = 0; epoch < epochs; ++epoch) {
        Rcpp::checkUserInterrupt();
        std::iota(order.begin(), order.end(), std::size_t(0));
        if (shuffle) {
            for (std::size_t remaining = count; remaining > 1; --remaining) {
                if ((remaining & 16383) == 0) Rcpp::checkUserInterrupt();
                std::swap(order[remaining - 1], order[bounded(rng, remaining)]);
            }
        }
        for (std::size_t pos = 0; pos < count; ++pos) {
            if ((pos & 16383) == 0) Rcpp::checkUserInterrupt();
            const Pair& p = pairs[order[pos]];
            // Test before multiplying so a very large positive product cannot overflow.
            const double mu = p.weight == 0 ? 0 : (rates[epoch] >= 1.0 / p.weight ? 1.0 : rates[epoch] * p.weight);
            const double mu_reverse = p.reverse == 0 ? 0 : (rates[epoch] >= 1.0 / p.reverse ? 1.0 : rates[epoch] * p.reverse);
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
                x(p.j, d) += (distance - p.target) * (0.5 * mu_reverse) * direction[d];
                if (!std::isfinite(x(p.i, d)) || !std::isfinite(x(p.j, d)))
                    Rcpp::stop("Nonfinite SGD coordinate update");
            }
        }
        updates += static_cast<double>(count); // Includes examined zero/zero pairs.
        completed = offset + epoch + 1;
        const bool scheduled=completed % checkpoint_every==0 || completed==schedule_epochs;
        if (scheduled || epoch+1==epochs) {
            const Score current=score(x,pairs,retain_best);
            const double diagnostic=retain_best ? current.profiled : current.raw;
            double motion=NA_REAL,residual=NA_REAL,relative=NA_REAL;
            if (scheduled && diagnostic<best_loss) {
                std::copy(x.begin(),x.end(),best.begin());best_loss=diagnostic;best_epoch=completed;
            }
            if (assess && scheduled) {
                motion=centered_motion(x,previous);residual=update_residual(x,pairs);
                relative=std::abs(diagnostic-previous_score)/std::max(std::abs(previous_score),1e-12);
                ++assessments;
                stable_count = completed>=min_epochs && relative<=rel_tol && motion<=motion_tol && residual<=residual_tol ? stable_count+1 : 0;
                previous=Rcpp::clone(x);previous_score=diagnostic;
                stabilized=stable_count>=patience;
            }
            record(completed,current,rates[epoch],scheduled,motion,residual,relative);
            if(stabilized) break;
        }
    }
    Rcpp::List next;
    if (keep_state) {
        std::ostringstream stream;stream<<rng;
        next=Rcpp::List::create(Rcpp::Named("coords")=Rcpp::clone(x),Rcpp::Named("epoch")=completed,
          Rcpp::Named("rng")=stream.str(),Rcpp::Named("best")=Rcpp::clone(best),
          Rcpp::Named("best_loss")=best_loss,Rcpp::Named("best_epoch")=best_epoch,
          Rcpp::Named("checkpoint_coords")=previous,Rcpp::Named("checkpoint_score")=previous_score,
          Rcpp::Named("stable_count")=stable_count,Rcpp::Named("assessments")=assessments);
    }
    // A chunk-boundary candidate may be returned, but must not affect the
    // canonical best checkpoint used when the same run is resumed later.
    const Score terminal=score(x,pairs,retain_best);
    const double terminal_loss=retain_best ? terminal.profiled : terminal.raw;
    const double diagnostic_best=std::min(best_loss,terminal_loss);
    const int diagnostic_epoch=terminal_loss<best_loss ? completed : best_epoch;
    if (retain_best && terminal_loss<best_loss) { best=Rcpp::clone(x);best_loss=terminal_loss;best_epoch=completed; }
    if (!retain_best) { best=x;best_epoch=completed;best_loss=NA_REAL; }
    return Rcpp::List::create(
        Rcpp::Named("conf")=best,Rcpp::Named("terminal_conf")=x,Rcpp::Named("niter")=completed,
        Rcpp::Named("stress")=best_loss,Rcpp::Named("best_epoch")=best_epoch,
        Rcpp::Named("pair_updates")=updates,Rcpp::Named("native_seed")=seed,
        Rcpp::Named("workspace_bytes")=static_cast<double>(workspace),
        Rcpp::Named("state")=next,Rcpp::Named("stabilized")=stabilized,
        Rcpp::Named("assessments")=assessments,Rcpp::Named("best_diagnostic")=diagnostic_best,
        Rcpp::Named("best_diagnostic_epoch")=diagnostic_epoch,
        Rcpp::Named("trace")=Rcpp::DataFrame::create(
          Rcpp::Named("epoch")=epoch_trace,Rcpp::Named("raw_stress")=raw_trace,
          Rcpp::Named("profiled_stress")=profiled_trace,Rcpp::Named("pair_updates")=updates_trace,
          Rcpp::Named("learning_rate")=rate_trace,Rcpp::Named("scheduled")=scheduled_trace,
          Rcpp::Named("relative_change")=relative_trace,Rcpp::Named("centered_motion")=motion_trace,
          Rcpp::Named("update_residual")=residual_trace,Rcpp::Named("stable_checkpoints")=streak_trace));
}

// Integrity check for serialized continuation data; not a security primitive.
// [[Rcpp::export]]
std::string grip_sgd_state_hash_cpp(Rcpp::RawVector bytes) {
    std::uint64_t value=14695981039346656037ULL;
    for(Rbyte byte:bytes) { value^=byte;value*=1099511628211ULL; }
    std::ostringstream out;out<<std::hex<<value;return out.str();
}
