#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct GramState {
    double energy;
    double edgeEnergy;
    double gramEnergy;
    double gradNorm2;
    double edgeScale;
    double edgeRelRmse;
    double gramRelRmse;
    Rcpp::NumericMatrix gradient;
};

void recenter_matrix(Rcpp::NumericMatrix &coords)
{
    int n = coords.nrow();
    int dim = coords.ncol();
    for(int j = 0; j < dim; j++){
        double mean = 0.0;
        for(int i = 0; i < n; i++)
            mean += coords(i, j);
        mean /= static_cast<double>(n);
        for(int i = 0; i < n; i++)
            coords(i, j) -= mean;
    }
}

void validate_common_args(int max_iter,
                          double lambda_edge,
                          double lambda_gram,
                          double edge_length_epsilon,
                          double initial_step,
                          double step_shrink,
                          double armijo_factor,
                          double grad_tol,
                          double min_step,
                          double distance_floor)
{
    if(max_iter < 0)
        Rcpp::stop("max_iter must be a non-negative integer");
    if(!std::isfinite(lambda_edge) || lambda_edge < 0.0)
        Rcpp::stop("lambda_edge must be finite and >= 0");
    if(!std::isfinite(lambda_gram) || lambda_gram < 0.0)
        Rcpp::stop("lambda_gram must be finite and >= 0");
    if(!std::isfinite(edge_length_epsilon) || edge_length_epsilon < 0.0)
        Rcpp::stop("edge_length_epsilon must be finite and >= 0");
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
    if(!std::isfinite(distance_floor) || distance_floor <= 0.0)
        Rcpp::stop("distance_floor must be finite and > 0");
}

void validate_inputs(const Rcpp::IntegerMatrix &edges,
                     const Rcpp::NumericVector &edge_weights,
                     const Rcpp::NumericVector &edge_stiffness,
                     const Rcpp::IntegerVector &star_center,
                     const Rcpp::IntegerVector &star_v1,
                     const Rcpp::IntegerVector &star_v2,
                     const Rcpp::NumericVector &star_w1,
                     const Rcpp::NumericVector &star_w2,
                     const Rcpp::NumericVector &star_cos,
                     const Rcpp::NumericVector &star_weight,
                     const Rcpp::NumericMatrix &coords)
{
    if(coords.ncol() != 2 && coords.ncol() != 3)
        Rcpp::stop("coords must have 2 or 3 columns");
    if(coords.nrow() < 2)
        Rcpp::stop("coords must have at least two rows");
    for(int i = 0; i < coords.size(); i++){
        if(!std::isfinite(coords[i]))
            Rcpp::stop("coords must contain only finite values");
    }
    if(edges.ncol() != 2)
        Rcpp::stop("edges must be a two-column matrix");
    if(edge_weights.size() != edges.nrow() || edge_stiffness.size() != edges.nrow())
        Rcpp::stop("edge_weights and edge_stiffness must be parallel to edges");
    for(int e = 0; e < edges.nrow(); e++){
        int u = edges(e, 0);
        int v = edges(e, 1);
        if(u <= 0 || u > coords.nrow() || v <= 0 || v > coords.nrow() || u == v)
            Rcpp::stop("edges must be valid 1-based non-self indices");
        if(!std::isfinite(edge_weights[e]) || edge_weights[e] <= 0.0)
            Rcpp::stop("edge_weights must contain finite values > 0");
        if(!std::isfinite(edge_stiffness[e]) || edge_stiffness[e] <= 0.0)
            Rcpp::stop("edge_stiffness must contain finite values > 0");
    }
    int np = star_center.size();
    if(star_v1.size() != np || star_v2.size() != np ||
       star_w1.size() != np || star_w2.size() != np ||
       star_cos.size() != np || star_weight.size() != np)
        Rcpp::stop("star vectors must have equal length");
    for(int r = 0; r < np; r++){
        int u = star_center[r];
        int v = star_v1[r];
        int w = star_v2[r];
        if(u <= 0 || u > coords.nrow() || v <= 0 || v > coords.nrow() || w <= 0 || w > coords.nrow())
            Rcpp::stop("star indices must be valid 1-based vertex ids");
        if(u == v || u == w || v == w)
            Rcpp::stop("star pairs must contain three distinct vertices");
        if(!std::isfinite(star_w1[r]) || star_w1[r] <= 0.0 ||
           !std::isfinite(star_w2[r]) || star_w2[r] <= 0.0)
            Rcpp::stop("star edge lengths must contain finite values > 0");
        if(!std::isfinite(star_cos[r]) || star_cos[r] < -1.0 || star_cos[r] > 1.0)
            Rcpp::stop("star cosines must be finite values in [-1, 1]");
        if(!std::isfinite(star_weight[r]) || star_weight[r] < 0.0)
            Rcpp::stop("star weights must contain finite values >= 0");
    }
}

double fit_edge_scale(const Rcpp::IntegerMatrix &edges,
                      const Rcpp::NumericVector &edge_weights,
                      const Rcpp::NumericVector &edge_stiffness,
                      const Rcpp::NumericMatrix &coords,
                      double edge_length_epsilon,
                      double distance_floor)
{
    double numerator = 0.0;
    double denominator = 0.0;
    double eps2 = edge_length_epsilon * edge_length_epsilon;
    for(int e = 0; e < edges.nrow(); e++){
        int u = edges(e, 0) - 1;
        int v = edges(e, 1) - 1;
        double d2 = eps2;
        for(int j = 0; j < coords.ncol(); j++){
            double diff = coords(u, j) - coords(v, j);
            d2 += diff * diff;
        }
        double observed = std::sqrt(d2);
        double target = edge_weights[e];
        double k = edge_stiffness[e];
        if(std::isfinite(observed) && std::isfinite(target) && target > distance_floor &&
           std::isfinite(k) && k > 0.0){
            numerator += k * observed * target;
            denominator += k * target * target;
        }
    }
    if(!std::isfinite(denominator) || denominator <= 0.0)
        return NA_REAL;
    return numerator / denominator;
}

GramState evaluate_state(const Rcpp::IntegerMatrix &edges,
                         const Rcpp::NumericVector &edge_weights,
                         const Rcpp::NumericVector &edge_stiffness,
                         const Rcpp::IntegerVector &star_center,
                         const Rcpp::IntegerVector &star_v1,
                         const Rcpp::IntegerVector &star_v2,
                         const Rcpp::NumericVector &star_w1,
                         const Rcpp::NumericVector &star_w2,
                         const Rcpp::NumericVector &star_cos,
                         const Rcpp::NumericVector &star_weight,
                         const Rcpp::NumericMatrix &coords,
                         const std::string &scale_mode,
                         double fixed_scale,
                         double lambda_edge,
                         double lambda_gram,
                         double edge_length_epsilon,
                         double distance_floor)
{
    int n = coords.nrow();
    int dim = coords.ncol();
    GramState state;
    state.energy = 0.0;
    state.edgeEnergy = 0.0;
    state.gramEnergy = 0.0;
    state.gradNorm2 = 0.0;
    state.edgeRelRmse = NA_REAL;
    state.gramRelRmse = NA_REAL;
    state.gradient = Rcpp::NumericMatrix(n, dim);

    double edgeScale = 1.0;
    if(scale_mode == "profiled"){
        edgeScale = fit_edge_scale(edges, edge_weights, edge_stiffness, coords, edge_length_epsilon, distance_floor);
        if(!std::isfinite(edgeScale))
            edgeScale = 1.0;
    } else if(scale_mode == "identity") {
        edgeScale = 1.0;
    } else if(scale_mode == "user") {
        if(!std::isfinite(fixed_scale) || fixed_scale <= 0.0)
            Rcpp::stop("scale must be finite and > 0");
        edgeScale = fixed_scale;
    } else {
        Rcpp::stop("unknown scale_mode");
    }
    state.edgeScale = edgeScale;

    double eps2 = edge_length_epsilon * edge_length_epsilon;
    double edgeDenom = 0.0;
    double edgeResid2 = 0.0;
    for(int e = 0; e < edges.nrow(); e++){
        int u = edges(e, 0) - 1;
        int v = edges(e, 1) - 1;
        double d2 = eps2;
        for(int j = 0; j < dim; j++){
            double diff = coords(u, j) - coords(v, j);
            d2 += diff * diff;
        }
        double len = std::sqrt(d2);
        double target = edgeScale * edge_weights[e];
        double residual = len - target;
        double k = lambda_edge * edge_stiffness[e];
        state.edgeEnergy += 0.5 * k * residual * residual;
        edgeResid2 += edge_stiffness[e] * residual * residual;
        edgeDenom += edge_stiffness[e] * target * target;
        if(k > 0.0){
            double coeff = k * residual / len;
            for(int j = 0; j < dim; j++){
                double contribution = coeff * (coords(u, j) - coords(v, j));
                state.gradient(u, j) += contribution;
                state.gradient(v, j) -= contribution;
            }
        }
    }
    if(std::isfinite(edgeDenom) && edgeDenom > 0.0)
        state.edgeRelRmse = std::sqrt(edgeResid2 / edgeDenom);

    double gramDenom = 0.0;
    double gramResid2 = 0.0;
    double scale2 = edgeScale * edgeScale;
    for(int r = 0; r < star_center.size(); r++){
        double baseWeight = star_weight[r];
        if(baseWeight <= 0.0 || lambda_gram <= 0.0)
            continue;
        int u = star_center[r] - 1;
        int v = star_v1[r] - 1;
        int w = star_v2[r] - 1;
        double observed = 0.0;
        for(int j = 0; j < dim; j++)
            observed += (coords(v, j) - coords(u, j)) * (coords(w, j) - coords(u, j));
        double target = scale2 * star_w1[r] * star_w2[r] * star_cos[r];
        double residual = observed - target;
        double k = lambda_gram * baseWeight;
        state.gramEnergy += 0.5 * k * residual * residual;
        gramResid2 += k * residual * residual;
        gramDenom += k * target * target;
        double coeff = k * residual;
        for(int j = 0; j < dim; j++){
            double a = coords(v, j) - coords(u, j);
            double b = coords(w, j) - coords(u, j);
            state.gradient(v, j) += coeff * b;
            state.gradient(w, j) += coeff * a;
            state.gradient(u, j) -= coeff * (a + b);
        }
    }
    if(std::isfinite(gramDenom) && gramDenom > 0.0)
        state.gramRelRmse = std::sqrt(gramResid2 / gramDenom);
    state.energy = state.edgeEnergy + state.gramEnergy;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < dim; j++)
            state.gradNorm2 += state.gradient(i, j) * state.gradient(i, j);
    }
    return state;
}

Rcpp::DataFrame build_trace_df(const std::vector<int> &iteration,
                               const std::vector<double> &energy,
                               const std::vector<double> &edge_energy,
                               const std::vector<double> &gram_energy,
                               const std::vector<double> &gradient_norm,
                               const std::vector<double> &step,
                               const std::vector<bool> &accepted,
                               const std::vector<double> &edge_scale,
                               const std::vector<double> &edge_rel_rmse,
                               const std::vector<double> &gram_rel_rmse)
{
    return Rcpp::DataFrame::create(
        Rcpp::_["iteration"] = iteration,
        Rcpp::_["energy"] = energy,
        Rcpp::_["edge.energy"] = edge_energy,
        Rcpp::_["gram.energy"] = gram_energy,
        Rcpp::_["gradient_norm"] = gradient_norm,
        Rcpp::_["step"] = step,
        Rcpp::_["accepted"] = accepted,
        Rcpp::_["edge.scale"] = edge_scale,
        Rcpp::_["edge.rel.rmse"] = edge_rel_rmse,
        Rcpp::_["gram.rel.rmse"] = gram_rel_rmse,
        Rcpp::_["stringsAsFactors"] = false
    );
}

} // namespace

// [[Rcpp::export]]
Rcpp::List grip_optimize_kernel_gram_gkk_layout_cpp(
    Rcpp::IntegerMatrix edges,
    Rcpp::NumericVector edge_weights,
    Rcpp::NumericVector edge_stiffness,
    Rcpp::IntegerVector star_center,
    Rcpp::IntegerVector star_v1,
    Rcpp::IntegerVector star_v2,
    Rcpp::NumericVector star_w1,
    Rcpp::NumericVector star_w2,
    Rcpp::NumericVector star_cos,
    Rcpp::NumericVector star_weight,
    Rcpp::NumericMatrix coords,
    int max_iter,
    std::string scale_mode,
    double scale,
    double lambda_edge,
    double lambda_gram,
    double edge_length_epsilon,
    double initial_step,
    double step_shrink,
    double armijo_factor,
    double grad_tol,
    double min_step,
    double distance_floor,
    bool recenter,
    bool return_trace)
{
    validate_common_args(max_iter,
                         lambda_edge,
                         lambda_gram,
                         edge_length_epsilon,
                         initial_step,
                         step_shrink,
                         armijo_factor,
                         grad_tol,
                         min_step,
                         distance_floor);
    validate_inputs(edges,
                    edge_weights,
                    edge_stiffness,
                    star_center,
                    star_v1,
                    star_v2,
                    star_w1,
                    star_w2,
                    star_cos,
                    star_weight,
                    coords);
    if(recenter)
        recenter_matrix(coords);

    std::vector<int> trace_iteration;
    std::vector<double> trace_energy;
    std::vector<double> trace_edge_energy;
    std::vector<double> trace_gram_energy;
    std::vector<double> trace_gradient_norm;
    std::vector<double> trace_step;
    std::vector<bool> trace_accepted;
    std::vector<double> trace_edge_scale;
    std::vector<double> trace_edge_rel_rmse;
    std::vector<double> trace_gram_rel_rmse;
    Rcpp::List frames;
    frames.push_back(Rcpp::clone(coords));

    GramState state = evaluate_state(edges,
                                     edge_weights,
                                     edge_stiffness,
                                     star_center,
                                     star_v1,
                                     star_v2,
                                     star_w1,
                                     star_w2,
                                     star_cos,
                                     star_weight,
                                     coords,
                                     scale_mode,
                                     scale,
                                     lambda_edge,
                                     lambda_gram,
                                     edge_length_epsilon,
                                     distance_floor);
    trace_iteration.push_back(0);
    trace_energy.push_back(state.energy);
    trace_edge_energy.push_back(state.edgeEnergy);
    trace_gram_energy.push_back(state.gramEnergy);
    trace_gradient_norm.push_back(std::sqrt(state.gradNorm2));
    trace_step.push_back(NA_REAL);
    trace_accepted.push_back(true);
    trace_edge_scale.push_back(state.edgeScale);
    trace_edge_rel_rmse.push_back(state.edgeRelRmse);
    trace_gram_rel_rmse.push_back(state.gramRelRmse);

    for(int iter = 1; iter <= max_iter; iter++){
        double gradNorm = std::sqrt(state.gradNorm2);
        if(!std::isfinite(gradNorm) || gradNorm <= grad_tol)
            break;
        double step = initial_step;
        bool accepted = false;
        Rcpp::NumericMatrix candidate = Rcpp::clone(coords);
        GramState candidateState = state;
        while(std::isfinite(step) && step >= min_step){
            Rcpp::NumericMatrix proposal = Rcpp::clone(coords);
            for(int i = 0; i < proposal.nrow(); i++){
                for(int j = 0; j < proposal.ncol(); j++)
                    proposal(i, j) -= step * state.gradient(i, j);
            }
            if(recenter)
                recenter_matrix(proposal);
            GramState proposalState = evaluate_state(edges,
                                                     edge_weights,
                                                     edge_stiffness,
                                                     star_center,
                                                     star_v1,
                                                     star_v2,
                                                     star_w1,
                                                     star_w2,
                                                     star_cos,
                                                     star_weight,
                                                     proposal,
                                                     scale_mode,
                                                     scale,
                                                     lambda_edge,
                                                     lambda_gram,
                                                     edge_length_epsilon,
                                                     distance_floor);
            double targetEnergy = state.energy - armijo_factor * step * state.gradNorm2;
            if(std::isfinite(proposalState.energy) && proposalState.energy <= targetEnergy){
                candidate = proposal;
                candidateState = proposalState;
                accepted = true;
                break;
            }
            step *= step_shrink;
        }
        trace_iteration.push_back(iter);
        trace_energy.push_back(accepted ? candidateState.energy : state.energy);
        trace_edge_energy.push_back(accepted ? candidateState.edgeEnergy : state.edgeEnergy);
        trace_gram_energy.push_back(accepted ? candidateState.gramEnergy : state.gramEnergy);
        trace_gradient_norm.push_back(accepted ? std::sqrt(candidateState.gradNorm2) : std::sqrt(state.gradNorm2));
        trace_step.push_back(accepted ? step : NA_REAL);
        trace_accepted.push_back(accepted);
        trace_edge_scale.push_back(accepted ? candidateState.edgeScale : state.edgeScale);
        trace_edge_rel_rmse.push_back(accepted ? candidateState.edgeRelRmse : state.edgeRelRmse);
        trace_gram_rel_rmse.push_back(accepted ? candidateState.gramRelRmse : state.gramRelRmse);
        if(!accepted)
            break;
        coords = candidate;
        state = candidateState;
        if(return_trace)
            frames.push_back(Rcpp::clone(coords));
    }

    Rcpp::DataFrame trace = build_trace_df(trace_iteration,
                                           trace_energy,
                                           trace_edge_energy,
                                           trace_gram_energy,
                                           trace_gradient_norm,
                                           trace_step,
                                           trace_accepted,
                                           trace_edge_scale,
                                           trace_edge_rel_rmse,
                                           trace_gram_rel_rmse);
    if(!return_trace){
        frames = Rcpp::List::create(Rcpp::clone(coords));
    }
    return Rcpp::List::create(
        Rcpp::_["coords"] = coords,
        Rcpp::_["trace"] = trace,
        Rcpp::_["frames"] = frames
    );
}
