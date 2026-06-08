#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct EdgeState {
    double energy;
    double gradNorm2;
    double edgeScale;
    double edgeRelRmse;
    Rcpp::NumericMatrix gradient;
};

void validate_args(int max_iter,
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
                     const Rcpp::NumericMatrix &stiffness_matrix,
                     const Rcpp::NumericVector &mix_schedule,
                     const Rcpp::NumericMatrix &coords)
{
    if(coords.ncol() < 2)
        Rcpp::stop("coords must have at least two columns");
    if(coords.nrow() < 2)
        Rcpp::stop("coords must have at least two rows");
    for(int i = 0; i < coords.size(); i++){
        if(!std::isfinite(coords[i]))
            Rcpp::stop("coords must contain only finite values");
    }
    if(edges.ncol() != 2)
        Rcpp::stop("edges must be a two-column matrix");
    if(edge_weights.size() != edges.nrow())
        Rcpp::stop("edge_weights must be parallel to edges");
    if(stiffness_matrix.nrow() != edges.nrow())
        Rcpp::stop("stiffness_matrix must have one row per edge");
    if(stiffness_matrix.ncol() != mix_schedule.size())
        Rcpp::stop("stiffness_matrix columns must match mix_schedule length");
    for(int e = 0; e < edges.nrow(); e++){
        int u = edges(e, 0);
        int v = edges(e, 1);
        if(u <= 0 || u > coords.nrow() || v <= 0 || v > coords.nrow())
            Rcpp::stop("edges must be 1-based and within [1, nrow(coords)]");
        if(u == v)
            Rcpp::stop("edges must not contain self loops");
        if(!std::isfinite(edge_weights[e]) || edge_weights[e] <= 0.0)
            Rcpp::stop("edge_weights must contain finite values > 0");
    }
    for(int j = 0; j < stiffness_matrix.ncol(); j++){
        if(!std::isfinite(mix_schedule[j]) || mix_schedule[j] < 0.0 || mix_schedule[j] > 1.0)
            Rcpp::stop("mix_schedule must contain values in [0, 1]");
        for(int e = 0; e < stiffness_matrix.nrow(); e++){
            if(!std::isfinite(stiffness_matrix(e, j)) || stiffness_matrix(e, j) <= 0.0)
                Rcpp::stop("stiffness_matrix must contain finite values > 0");
        }
    }
}

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

double fit_edge_scale(const Rcpp::IntegerMatrix &edges,
                      const Rcpp::NumericVector &edge_weights,
                      const Rcpp::NumericVector &stiffness,
                      const Rcpp::NumericMatrix &coords,
                      double edge_length_epsilon,
                      double distance_floor)
{
    double numerator = 0.0;
    double denominator = 0.0;
    double eps2 = edge_length_epsilon * edge_length_epsilon;
    for(int e = 0; e < edges.nrow(); e++){
        double target = edge_weights[e];
        double k = stiffness[e];
        if(!std::isfinite(target) || target <= distance_floor || !std::isfinite(k) || k <= 0.0)
            continue;
        int u = edges(e, 0) - 1;
        int v = edges(e, 1) - 1;
        double d2 = eps2;
        for(int j = 0; j < coords.ncol(); j++){
            double diff = coords(u, j) - coords(v, j);
            d2 += diff * diff;
        }
        double observed = std::sqrt(d2);
        if(!std::isfinite(observed))
            continue;
        numerator += k * observed * target;
        denominator += k * target * target;
    }
    if(!std::isfinite(denominator) || denominator <= 0.0)
        return NA_REAL;
    return numerator / denominator;
}

EdgeState evaluate_state(const Rcpp::IntegerMatrix &edges,
                         const Rcpp::NumericVector &edge_weights,
                         const Rcpp::NumericVector &stiffness,
                         const Rcpp::NumericMatrix &coords,
                         const std::string &scale_mode,
                         double fixed_scale,
                         double edge_length_epsilon,
                         double distance_floor)
{
    int n = coords.nrow();
    int dim = coords.ncol();
    EdgeState state;
    state.energy = 0.0;
    state.gradNorm2 = 0.0;
    state.edgeRelRmse = NA_REAL;
    state.gradient = Rcpp::NumericMatrix(n, dim);

    double edgeScale = 1.0;
    if(scale_mode == "profiled"){
        edgeScale = fit_edge_scale(edges,
                                   edge_weights,
                                   stiffness,
                                   coords,
                                   edge_length_epsilon,
                                   distance_floor);
        if(!std::isfinite(edgeScale))
            edgeScale = 1.0;
    } else if(scale_mode == "identity") {
        edgeScale = 1.0;
    } else if(scale_mode == "fixed" || scale_mode == "user") {
        if(!std::isfinite(fixed_scale) || fixed_scale <= 0.0)
            Rcpp::stop("scale must be finite and > 0");
        edgeScale = fixed_scale;
    } else {
        Rcpp::stop("unknown scale_mode");
    }
    state.edgeScale = edgeScale;

    double eps2 = edge_length_epsilon * edge_length_epsilon;
    double denom = 0.0;
    double resid2 = 0.0;
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
        double k = stiffness[e];
        state.energy += 0.5 * k * residual * residual;
        resid2 += k * residual * residual;
        denom += k * target * target;
        double coeff = k * residual / len;
        for(int j = 0; j < dim; j++){
            double contribution = coeff * (coords(u, j) - coords(v, j));
            state.gradient(u, j) += contribution;
            state.gradient(v, j) -= contribution;
        }
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < dim; j++)
            state.gradNorm2 += state.gradient(i, j) * state.gradient(i, j);
    }
    if(std::isfinite(denom) && denom > 0.0)
        state.edgeRelRmse = std::sqrt(resid2 / denom);
    return state;
}

Rcpp::DataFrame build_trace_df(const std::vector<int> &stage,
                               const std::vector<double> &mix,
                               const std::vector<int> &iteration,
                               const std::vector<double> &energy,
                               const std::vector<double> &gradient_norm,
                               const std::vector<double> &step,
                               const std::vector<bool> &accepted,
                               const std::vector<double> &edge_scale,
                               const std::vector<double> &edge_rel_rmse)
{
    return Rcpp::DataFrame::create(
        Rcpp::_["stage"] = stage,
        Rcpp::_["mix"] = mix,
        Rcpp::_["iteration"] = iteration,
        Rcpp::_["energy"] = energy,
        Rcpp::_["gradient_norm"] = gradient_norm,
        Rcpp::_["step"] = step,
        Rcpp::_["accepted"] = accepted,
        Rcpp::_["edge.scale"] = edge_scale,
        Rcpp::_["edge.rel.rmse"] = edge_rel_rmse,
        Rcpp::_["stringsAsFactors"] = false
    );
}

} // namespace

// [[Rcpp::export]]
Rcpp::List grip_optimize_edge_isometric_layout_cpp(
    Rcpp::IntegerMatrix edges,
    Rcpp::NumericVector edge_weights,
    Rcpp::NumericMatrix stiffness_matrix,
    Rcpp::NumericVector mix_schedule,
    Rcpp::NumericMatrix coords,
    int max_iter,
    std::string scale_mode,
    double scale,
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
    validate_args(max_iter,
                  edge_length_epsilon,
                  initial_step,
                  step_shrink,
                  armijo_factor,
                  grad_tol,
                  min_step,
                  distance_floor);
    validate_inputs(edges, edge_weights, stiffness_matrix, mix_schedule, coords);
    if(recenter)
        recenter_matrix(coords);

    std::vector<int> trace_stage;
    std::vector<double> trace_mix;
    std::vector<int> trace_iteration;
    std::vector<double> trace_energy;
    std::vector<double> trace_gradient_norm;
    std::vector<double> trace_step;
    std::vector<bool> trace_accepted;
    std::vector<double> trace_edge_scale;
    std::vector<double> trace_edge_rel_rmse;
    std::vector<int> final_stage;
    std::vector<double> final_mix;
    std::vector<int> final_iteration;
    std::vector<double> final_energy;
    std::vector<double> final_gradient_norm;
    std::vector<double> final_step;
    std::vector<bool> final_accepted;
    std::vector<double> final_edge_scale;
    std::vector<double> final_edge_rel_rmse;
    Rcpp::List frames;
    if(return_trace)
        frames.push_back(Rcpp::clone(coords));

    double gradTol2 = grad_tol * grad_tol;

    for(int stage = 0; stage < stiffness_matrix.ncol(); stage++){
        Rcpp::NumericVector stiffness = stiffness_matrix(Rcpp::_, stage);
        EdgeState state = evaluate_state(edges,
                                         edge_weights,
                                         stiffness,
                                         coords,
                                         scale_mode,
                                         scale,
                                         edge_length_epsilon,
                                         distance_floor);
        int currentIteration = 0;
        double currentStep = NA_REAL;
        bool currentAccepted = true;
        if(return_trace){
            trace_stage.push_back(stage + 1);
            trace_mix.push_back(mix_schedule[stage]);
            trace_iteration.push_back(currentIteration);
            trace_energy.push_back(state.energy);
            trace_gradient_norm.push_back(std::sqrt(state.gradNorm2));
            trace_step.push_back(currentStep);
            trace_accepted.push_back(currentAccepted);
            trace_edge_scale.push_back(state.edgeScale);
            trace_edge_rel_rmse.push_back(state.edgeRelRmse);
        }

        for(int iter = 1; iter <= max_iter; iter++){
            if(!std::isfinite(state.energy) || state.gradNorm2 <= gradTol2)
                break;

            double step = initial_step;
            bool accepted = false;
            Rcpp::NumericMatrix candidate = Rcpp::clone(coords);
            EdgeState candidateState = state;

            while(std::isfinite(step) && step >= min_step){
                Rcpp::NumericMatrix proposal = Rcpp::clone(coords);
                for(int i = 0; i < proposal.nrow(); i++){
                    for(int j = 0; j < proposal.ncol(); j++)
                        proposal(i, j) -= step * state.gradient(i, j);
                }
                if(recenter)
                    recenter_matrix(proposal);

                EdgeState proposalState = evaluate_state(edges,
                                                         edge_weights,
                                                         stiffness,
                                                         proposal,
                                                         scale_mode,
                                                         scale,
                                                         edge_length_epsilon,
                                                         distance_floor);
                double targetEnergy = state.energy - armijo_factor * step * state.gradNorm2;
                if(std::isfinite(proposalState.energy) && proposalState.energy <= targetEnergy){
                    accepted = true;
                    candidate = proposal;
                    candidateState = proposalState;
                    break;
                }
                step *= step_shrink;
            }

            currentIteration = iter;
            currentStep = accepted ? step : NA_REAL;
            currentAccepted = accepted;
            if(return_trace){
                trace_stage.push_back(stage + 1);
                trace_mix.push_back(mix_schedule[stage]);
                trace_iteration.push_back(currentIteration);
                trace_energy.push_back(accepted ? candidateState.energy : state.energy);
                trace_gradient_norm.push_back(std::sqrt(accepted ? candidateState.gradNorm2 : state.gradNorm2));
                trace_step.push_back(currentStep);
                trace_accepted.push_back(currentAccepted);
                trace_edge_scale.push_back(accepted ? candidateState.edgeScale : state.edgeScale);
                trace_edge_rel_rmse.push_back(accepted ? candidateState.edgeRelRmse : state.edgeRelRmse);
            }

            if(!accepted)
                break;

            coords = candidate;
            state = candidateState;
            if(return_trace)
                frames.push_back(Rcpp::clone(coords));
        }

        final_stage.push_back(stage + 1);
        final_mix.push_back(mix_schedule[stage]);
        final_iteration.push_back(currentIteration);
        final_energy.push_back(state.energy);
        final_gradient_norm.push_back(std::sqrt(state.gradNorm2));
        final_step.push_back(currentStep);
        final_accepted.push_back(currentAccepted);
        final_edge_scale.push_back(state.edgeScale);
        final_edge_rel_rmse.push_back(state.edgeRelRmse);
    }

    Rcpp::RObject trace;
    Rcpp::RObject frame_output;
    if(return_trace){
        trace = build_trace_df(trace_stage,
                               trace_mix,
                               trace_iteration,
                               trace_energy,
                               trace_gradient_norm,
                               trace_step,
                               trace_accepted,
                               trace_edge_scale,
                               trace_edge_rel_rmse);
        frame_output = frames;
    } else {
        trace = R_NilValue;
        frame_output = R_NilValue;
    }

    return Rcpp::List::create(
        Rcpp::_["coords"] = coords,
        Rcpp::_["trace"] = trace,
        Rcpp::_["stage_final"] = build_trace_df(final_stage,
                                                final_mix,
                                                final_iteration,
                                                final_energy,
                                                final_gradient_norm,
                                                final_step,
                                                final_accepted,
                                                final_edge_scale,
                                                final_edge_rel_rmse),
        Rcpp::_["frames"] = frame_output
    );
}
