#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

struct SingleVertexState {
    double objective;
    double gradNorm2;
    std::vector<double> gradient;
};

void validate_insert_args(const Rcpp::NumericMatrix &anchor_coords,
                          const Rcpp::NumericVector &anchor_distance,
                          int max_iter,
                          double initial_step,
                          double step_shrink,
                          double armijo_factor,
                          double grad_tol,
                          double min_step)
{
    if(anchor_coords.nrow() <= 0)
        Rcpp::stop("anchor_coords must have at least one row");
    if(anchor_coords.ncol() <= 0)
        Rcpp::stop("anchor_coords must have at least one column");
    if(anchor_distance.size() != anchor_coords.nrow())
        Rcpp::stop("anchor_distance length must match nrow(anchor_coords)");
    for(int i = 0; i < anchor_coords.nrow(); i++){
        for(int j = 0; j < anchor_coords.ncol(); j++){
            double value = anchor_coords(i, j);
            if(!std::isfinite(value))
                Rcpp::stop("anchor_coords must contain finite values");
        }
        double dist = anchor_distance[i];
        if(!std::isfinite(dist) || dist < 0.0)
            Rcpp::stop("anchor_distance must contain finite values >= 0");
    }
    if(max_iter < 0)
        Rcpp::stop("max_iter must be a non-negative integer");
    if(!std::isfinite(initial_step) || initial_step <= 0.0)
        Rcpp::stop("initial_step must be finite and > 0");
    if(!std::isfinite(step_shrink) || step_shrink <= 0.0 || step_shrink >= 1.0)
        Rcpp::stop("step_shrink must be finite and in (0, 1)");
    if(!std::isfinite(armijo_factor) || armijo_factor < 0.0)
        Rcpp::stop("armijo_factor must be finite and >= 0");
    if(!std::isfinite(grad_tol) || grad_tol < 0.0)
        Rcpp::stop("grad_tol must be finite and >= 0");
    if(!std::isfinite(min_step) || min_step <= 0.0)
        Rcpp::stop("min_step must be finite and > 0");
}

std::vector<double> validate_weights(const Rcpp::Nullable<Rcpp::NumericVector> &anchor_weights,
                                     int anchor_count)
{
    std::vector<double> weights(static_cast<size_t>(anchor_count), 1.0);
    if(anchor_weights.isNull())
        return weights;

    Rcpp::NumericVector raw = anchor_weights.get();
    if(raw.size() != anchor_count)
        Rcpp::stop("anchor_weights length must match nrow(anchor_coords)");
    for(int i = 0; i < anchor_count; i++){
        double value = raw[i];
        if(!std::isfinite(value) || value <= 0.0)
            Rcpp::stop("anchor_weights must contain finite values > 0");
        weights[static_cast<size_t>(i)] = value;
    }
    return weights;
}

bool solve_linear_system(std::vector<std::vector<double>> A,
                         std::vector<double> b,
                         std::vector<double> &x)
{
    const size_t n = A.size();
    if(n == 0 || b.size() != n)
        return false;
    for(size_t i = 0; i < n; i++){
        if(A[i].size() != n)
            return false;
    }

    for(size_t col = 0; col < n; col++){
        size_t pivot = col;
        double best = std::abs(A[pivot][col]);
        for(size_t row = col + 1; row < n; row++){
            double val = std::abs(A[row][col]);
            if(val > best){
                best = val;
                pivot = row;
            }
        }
        if(best <= 1e-12)
            return false;
        if(pivot != col){
            std::swap(A[pivot], A[col]);
            std::swap(b[pivot], b[col]);
        }
        const double diag = A[col][col];
        for(size_t row = col + 1; row < n; row++){
            const double factor = A[row][col] / diag;
            if(std::abs(factor) <= 0.0)
                continue;
            for(size_t k = col; k < n; k++)
                A[row][k] -= factor * A[col][k];
            b[row] -= factor * b[col];
        }
    }

    x.assign(n, 0.0);
    for(int row = static_cast<int>(n) - 1; row >= 0; row--){
        double rhs = b[static_cast<size_t>(row)];
        for(size_t col = static_cast<size_t>(row) + 1; col < n; col++)
            rhs -= A[static_cast<size_t>(row)][col] * x[col];
        const double diag = A[static_cast<size_t>(row)][static_cast<size_t>(row)];
        if(std::abs(diag) <= 1e-12)
            return false;
        x[static_cast<size_t>(row)] = rhs / diag;
    }
    return true;
}

std::vector<double> weighted_barycenter(const Rcpp::NumericMatrix &anchor_coords,
                                        const std::vector<double> &weights,
                                        const Rcpp::NumericVector &anchor_distance)
{
    const int anchor_count = anchor_coords.nrow();
    const int dim = anchor_coords.ncol();
    std::vector<double> center(static_cast<size_t>(dim), 0.0);
    double weight_sum = 0.0;
    for(int i = 0; i < anchor_count; i++){
        const double wi = weights[static_cast<size_t>(i)];
        weight_sum += wi;
        for(int j = 0; j < dim; j++)
            center[static_cast<size_t>(j)] += wi * anchor_coords(i, j);
    }
    if(weight_sum <= 0.0)
        weight_sum = 1.0;
    for(int j = 0; j < dim; j++)
        center[static_cast<size_t>(j)] /= weight_sum;

    if(anchor_count == 1 && dim >= 1)
        center[0] += anchor_distance[0];
    return center;
}

std::vector<double> linearized_trilateration_start(const Rcpp::NumericMatrix &anchor_coords,
                                                   const Rcpp::NumericVector &anchor_distance,
                                                   const std::vector<double> &weights)
{
    const int anchor_count = anchor_coords.nrow();
    const int dim = anchor_coords.ncol();
    std::vector<double> fallback = weighted_barycenter(anchor_coords, weights, anchor_distance);
    if(anchor_count < 2)
        return fallback;

    std::vector<std::vector<double>> normal(static_cast<size_t>(dim), std::vector<double>(static_cast<size_t>(dim), 0.0));
    std::vector<double> rhs(static_cast<size_t>(dim), 0.0);
    const int ref = 0;

    double ref_norm2 = 0.0;
    for(int j = 0; j < dim; j++)
        ref_norm2 += anchor_coords(ref, j) * anchor_coords(ref, j);
    const double ref_dist2 = anchor_distance[ref] * anchor_distance[ref];

    for(int i = 1; i < anchor_count; i++){
        std::vector<double> row(static_cast<size_t>(dim), 0.0);
        double ai_norm2 = 0.0;
        for(int j = 0; j < dim; j++){
            const double diff = 2.0 * (anchor_coords(i, j) - anchor_coords(ref, j));
            row[static_cast<size_t>(j)] = diff;
            ai_norm2 += anchor_coords(i, j) * anchor_coords(i, j);
        }
        const double dist2 = anchor_distance[i] * anchor_distance[i];
        const double b = ai_norm2 - ref_norm2 - dist2 + ref_dist2;
        const double wi = std::sqrt(std::max(weights[static_cast<size_t>(i)], 1e-12));

        for(int a = 0; a < dim; a++){
            rhs[static_cast<size_t>(a)] += wi * row[static_cast<size_t>(a)] * (wi * b);
            for(int c = 0; c < dim; c++)
                normal[static_cast<size_t>(a)][static_cast<size_t>(c)] += wi * row[static_cast<size_t>(a)] * (wi * row[static_cast<size_t>(c)]);
        }
    }

    for(int j = 0; j < dim; j++)
        normal[static_cast<size_t>(j)][static_cast<size_t>(j)] += 1e-8;

    std::vector<double> solution;
    if(!solve_linear_system(normal, rhs, solution))
        return fallback;
    return solution;
}

SingleVertexState evaluate_single_vertex(const std::vector<double> &x,
                                         const Rcpp::NumericMatrix &anchor_coords,
                                         const Rcpp::NumericVector &anchor_distance,
                                         const std::vector<double> &weights)
{
    const int anchor_count = anchor_coords.nrow();
    const int dim = anchor_coords.ncol();
    SingleVertexState state;
    state.objective = 0.0;
    state.gradNorm2 = 0.0;
    state.gradient.assign(static_cast<size_t>(dim), 0.0);

    for(int i = 0; i < anchor_count; i++){
        double norm2 = 0.0;
        std::vector<double> diff(static_cast<size_t>(dim), 0.0);
        for(int j = 0; j < dim; j++){
            diff[static_cast<size_t>(j)] = x[static_cast<size_t>(j)] - anchor_coords(i, j);
            norm2 += diff[static_cast<size_t>(j)] * diff[static_cast<size_t>(j)];
        }
        const double radius = std::sqrt(std::max(norm2, 1e-16));
        const double residual = radius - anchor_distance[i];
        const double wi = weights[static_cast<size_t>(i)];
        state.objective += 0.5 * wi * residual * residual;
        const double scale = wi * residual / radius;
        for(int j = 0; j < dim; j++)
            state.gradient[static_cast<size_t>(j)] += scale * diff[static_cast<size_t>(j)];
    }

    for(size_t j = 0; j < state.gradient.size(); j++)
        state.gradNorm2 += state.gradient[j] * state.gradient[j];
    return state;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List grip_geodesic_misf_insert_vertex_cpp(Rcpp::NumericMatrix anchor_coords,
                                                Rcpp::NumericVector anchor_distance,
                                                Rcpp::Nullable<Rcpp::NumericVector> anchor_weights = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericVector> init_coord = R_NilValue,
                                                int max_iter = 64,
                                                double initial_step = 1.0,
                                                double step_shrink = 0.5,
                                                double armijo_factor = 1e-4,
                                                double grad_tol = 1e-8,
                                                double min_step = 1e-8)
{
    validate_insert_args(
        anchor_coords,
        anchor_distance,
        max_iter,
        initial_step,
        step_shrink,
        armijo_factor,
        grad_tol,
        min_step
    );

    const int anchor_count = anchor_coords.nrow();
    const int dim = anchor_coords.ncol();
    std::vector<double> weights = validate_weights(anchor_weights, anchor_count);

    std::vector<double> current;
    if(init_coord.isNotNull()){
        Rcpp::NumericVector init = init_coord.get();
        if(init.size() != dim)
            Rcpp::stop("init_coord length must match ncol(anchor_coords)");
        current.assign(static_cast<size_t>(dim), 0.0);
        for(int j = 0; j < dim; j++){
            double value = init[j];
            if(!std::isfinite(value))
                Rcpp::stop("init_coord must contain finite values");
            current[static_cast<size_t>(j)] = value;
        }
    } else {
        current = linearized_trilateration_start(anchor_coords, anchor_distance, weights);
    }

    Rcpp::NumericVector initial_coord(dim);
    for(int j = 0; j < dim; j++)
        initial_coord[j] = current[static_cast<size_t>(j)];

    SingleVertexState state = evaluate_single_vertex(current, anchor_coords, anchor_distance, weights);
    const double initial_objective = state.objective;
    int iterations = 0;
    bool converged = (std::sqrt(state.gradNorm2) <= grad_tol);
    double last_step = 0.0;

    while(iterations < max_iter && !converged){
        const double grad_norm2 = state.gradNorm2;
        if(grad_norm2 <= grad_tol * grad_tol){
            converged = true;
            break;
        }

        double step = initial_step;
        bool accepted = false;
        while(step >= min_step){
            std::vector<double> trial = current;
            for(int j = 0; j < dim; j++)
                trial[static_cast<size_t>(j)] -= step * state.gradient[static_cast<size_t>(j)];
            SingleVertexState trial_state = evaluate_single_vertex(trial, anchor_coords, anchor_distance, weights);
            if(trial_state.objective <= state.objective - armijo_factor * step * grad_norm2){
                current.swap(trial);
                state = trial_state;
                accepted = true;
                last_step = step;
                break;
            }
            step *= step_shrink;
        }

        iterations += 1;
        if(!accepted)
            break;
        if(std::sqrt(state.gradNorm2) <= grad_tol)
            converged = true;
    }

    Rcpp::NumericVector coord(dim);
    for(int j = 0; j < dim; j++)
        coord[j] = current[static_cast<size_t>(j)];
    return Rcpp::List::create(
        Rcpp::Named("coord") = coord,
        Rcpp::Named("initial_coord") = initial_coord,
        Rcpp::Named("objective") = state.objective,
        Rcpp::Named("initial_objective") = initial_objective,
        Rcpp::Named("grad_norm") = std::sqrt(state.gradNorm2),
        Rcpp::Named("iterations") = iterations,
        Rcpp::Named("converged") = converged,
        Rcpp::Named("accepted_step") = last_step,
        Rcpp::Named("anchor_count") = anchor_count,
        Rcpp::Named("dim") = dim
    );
}
