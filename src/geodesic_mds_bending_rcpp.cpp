#include <Rcpp.h>
#include <cmath>
#include <vector>

#include "Point.h"

namespace {

struct FlatCache {
    std::vector<int> offsets, u, v;
    std::vector<double> coeff, graph;
};

struct BendState {
    double energy, gmdsEnergy, anchorEnergy, bendEnergy, gradNorm2;
    std::vector<Point<>> gradient;
};

void validate_args(int max_iter, double eps, double initial_step, double step_shrink,
                   double armijo_factor, double grad_tol, double min_step)
{
    if(max_iter < 0) Rcpp::stop("max_iter must be non-negative");
    if(!std::isfinite(eps) || eps < 0.0) Rcpp::stop("edge_length_epsilon must be finite and >= 0");
    if(!std::isfinite(initial_step) || initial_step <= 0.0) Rcpp::stop("initial_step must be finite and > 0");
    if(!std::isfinite(step_shrink) || step_shrink <= 0.0 || step_shrink >= 1.0) Rcpp::stop("step_shrink must be in (0,1)");
    if(!std::isfinite(armijo_factor) || armijo_factor < 0.0) Rcpp::stop("armijo_factor must be >= 0");
    if(!std::isfinite(grad_tol) || grad_tol < 0.0) Rcpp::stop("grad_tol must be >= 0");
    if(!std::isfinite(min_step) || min_step <= 0.0) Rcpp::stop("min_step must be > 0");
}

std::vector<Point<>> matrix_to_points(const Rcpp::NumericMatrix &coords)
{
    std::vector<Point<>> out(static_cast<size_t>(coords.nrow()));
    for(int i = 0; i < coords.nrow(); i++){
        out[static_cast<size_t>(i)].set_to_zero();
        out[static_cast<size_t>(i)].setX(coords(i, 0));
        if(coords.ncol() > 1) out[static_cast<size_t>(i)].setY(coords(i, 1));
        if(coords.ncol() > 2) out[static_cast<size_t>(i)].setZ(coords(i, 2));
    }
    return out;
}

Rcpp::NumericMatrix points_to_matrix(const std::vector<Point<>> &coords, int dim)
{
    Rcpp::NumericMatrix out(static_cast<int>(coords.size()), dim);
    for(int i = 0; i < out.nrow(); i++){
        out(i, 0) = coords[static_cast<size_t>(i)].getX();
        if(dim > 1) out(i, 1) = coords[static_cast<size_t>(i)].getY();
        if(dim > 2) out(i, 2) = coords[static_cast<size_t>(i)].getZ();
    }
    return out;
}

void recenter_points(std::vector<Point<>> &coords)
{
    if(coords.empty()) return;
    Point<> mean; mean.set_to_zero();
    for(size_t i = 0; i < coords.size(); i++) mean += coords[i];
    mean /= static_cast<double>(coords.size());
    for(size_t i = 0; i < coords.size(); i++) coords[i] -= mean;
}

std::vector<double> resolve_schedule(Rcpp::Nullable<Rcpp::NumericVector> weights, int max_iter)
{
    std::vector<double> out(static_cast<size_t>(max_iter + 1), 0.0);
    if(weights.isNull()) return out;
    Rcpp::NumericVector raw = weights.get();
    if(raw.size() == 0) return out;
    if(raw.size() != 1 && raw.size() != max_iter + 1) Rcpp::stop("schedule must have length 1 or max_iter + 1");
    for(int i = 0; i <= max_iter; i++){
        double value = raw.size() == 1 ? raw[0] : raw[i];
        if(!std::isfinite(value) || value < 0.0) Rcpp::stop("schedule entries must be finite and >= 0");
        out[static_cast<size_t>(i)] = value;
    }
    return out;
}

BendState evaluate_state(const std::vector<Point<>> &coords, const FlatCache &cache,
                         double eps2, const std::vector<Point<>> *anchor, double anchorWeight,
                         const std::vector<int> &bendA, const std::vector<int> &bendB,
                         const std::vector<int> &bendC, double bendWeight)
{
    BendState state{0.0, 0.0, 0.0, 0.0, 0.0, std::vector<Point<>>(coords.size())};
    for(size_t i = 0; i < state.gradient.size(); i++) state.gradient[i].set_to_zero();

    for(size_t pairIndex = 0; pairIndex < cache.graph.size(); pairIndex++){
        int begin = cache.offsets[pairIndex], end = cache.offsets[pairIndex + 1];
        double h = 0.0;
        for(int edgeIndex = begin; edgeIndex < end; edgeIndex++){
            Point<> diff = coords[static_cast<size_t>(cache.u[edgeIndex])] - coords[static_cast<size_t>(cache.v[edgeIndex])];
            h += cache.coeff[edgeIndex] * std::sqrt(diff.fnorm2() + eps2);
        }
        double resid = h - cache.graph[pairIndex];
        state.gmdsEnergy += 0.5 * resid * resid;
        for(int edgeIndex = begin; edgeIndex < end; edgeIndex++){
            Point<> diff = coords[static_cast<size_t>(cache.u[edgeIndex])] - coords[static_cast<size_t>(cache.v[edgeIndex])];
            double len = std::sqrt(diff.fnorm2() + eps2);
            if(len <= 0.0) continue;
            Point<> step = diff * (resid * cache.coeff[edgeIndex] / len);
            state.gradient[static_cast<size_t>(cache.u[edgeIndex])] += step;
            state.gradient[static_cast<size_t>(cache.v[edgeIndex])] -= step;
        }
    }

    if(anchor && anchorWeight > 0.0){
        double raw = 0.0;
        for(size_t i = 0; i < coords.size(); i++){
            Point<> diff = coords[i] - (*anchor)[i];
            raw += diff.fnorm2();
            state.gradient[i] += diff * (2.0 * anchorWeight);
        }
        state.anchorEnergy = anchorWeight * raw;
    }

    if(!bendA.empty() && bendWeight > 0.0){
        double raw = 0.0, scale = 2.0 * bendWeight / static_cast<double>(bendA.size());
        for(size_t i = 0; i < bendA.size(); i++){
            Point<> resid = coords[static_cast<size_t>(bendA[i])] - 2.0 * coords[static_cast<size_t>(bendB[i])] + coords[static_cast<size_t>(bendC[i])];
            raw += resid.fnorm2();
            state.gradient[static_cast<size_t>(bendA[i])] += resid * scale;
            state.gradient[static_cast<size_t>(bendB[i])] -= resid * (2.0 * scale);
            state.gradient[static_cast<size_t>(bendC[i])] += resid * scale;
        }
        state.bendEnergy = bendWeight * (raw / static_cast<double>(bendA.size()));
    }

    state.energy = state.gmdsEnergy + state.anchorEnergy + state.bendEnergy;
    for(size_t i = 0; i < state.gradient.size(); i++) state.gradNorm2 += state.gradient[i].fnorm2();
    return state;
}

Rcpp::DataFrame build_trace_df(const std::vector<int> &iteration, const std::vector<double> &energy,
                               const std::vector<double> &gmds_energy, const std::vector<double> &anchor_energy,
                               const std::vector<double> &bend_energy, const std::vector<double> &gradient_norm,
                               const std::vector<double> &step, const std::vector<bool> &accepted,
                               const std::vector<double> &anchor_weight, const std::vector<double> &bend_weight)
{
    return Rcpp::DataFrame::create(Rcpp::_["iteration"] = iteration, Rcpp::_["energy"] = energy,
                                   Rcpp::_["gmds_energy"] = gmds_energy, Rcpp::_["anchor_energy"] = anchor_energy,
                                   Rcpp::_["smooth_energy"] = Rcpp::NumericVector(iteration.size(), 0.0),
                                   Rcpp::_["bend_energy"] = bend_energy, Rcpp::_["gradient_norm"] = gradient_norm,
                                   Rcpp::_["step"] = step, Rcpp::_["accepted"] = accepted,
                                   Rcpp::_["anchor_weight"] = anchor_weight,
                                   Rcpp::_["smooth_weight"] = Rcpp::NumericVector(iteration.size(), 0.0),
                                   Rcpp::_["bend_weight"] = bend_weight, Rcpp::_["stringsAsFactors"] = false);
}

} // namespace

// [[Rcpp::export]]
Rcpp::List grip_optimize_geodesic_mds_flat_bending_cpp(
    Rcpp::IntegerVector flat_pair_edge_offsets, Rcpp::IntegerVector flat_edge_u,
    Rcpp::IntegerVector flat_edge_v, Rcpp::NumericVector flat_edge_coeff,
    Rcpp::NumericVector pair_graph_distance, Rcpp::NumericMatrix coords,
    int max_iter, double edge_length_epsilon, double initial_step, double step_shrink,
    double armijo_factor, double grad_tol, double min_step, bool recenter, bool return_trace,
    Rcpp::Nullable<Rcpp::NumericMatrix> anchor_coords = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> anchor_weights = R_NilValue,
    Rcpp::IntegerVector bend_a = Rcpp::IntegerVector(), Rcpp::IntegerVector bend_b = Rcpp::IntegerVector(),
    Rcpp::IntegerVector bend_c = Rcpp::IntegerVector(), Rcpp::Nullable<Rcpp::NumericVector> bend_weights = R_NilValue)
{
    if(coords.ncol() != 2 && coords.ncol() != 3) Rcpp::stop("coords must have 2 or 3 columns");
    validate_args(max_iter, edge_length_epsilon, initial_step, step_shrink, armijo_factor, grad_tol, min_step);
    if(flat_pair_edge_offsets.size() != pair_graph_distance.size() + 1) Rcpp::stop("flat_pair_edge_offsets must have length length(pair_graph_distance) + 1");
    if(flat_edge_u.size() != flat_edge_v.size() || flat_edge_u.size() != flat_edge_coeff.size()) Rcpp::stop("flat edge arrays must have the same length");
    if(bend_a.size() != bend_b.size() || bend_a.size() != bend_c.size()) Rcpp::stop("bend arrays must have the same length");

    FlatCache cache{Rcpp::as<std::vector<int>>(flat_pair_edge_offsets), Rcpp::as<std::vector<int>>(flat_edge_u),
                    Rcpp::as<std::vector<int>>(flat_edge_v), Rcpp::as<std::vector<double>>(flat_edge_coeff),
                    Rcpp::as<std::vector<double>>(pair_graph_distance)};
    std::vector<Point<>> current = matrix_to_points(coords), anchor;
    bool useAnchor = anchor_coords.isNotNull();
    if(useAnchor){ Rcpp::NumericMatrix a = anchor_coords.get(); if(a.nrow() != coords.nrow() || a.ncol() != coords.ncol()) Rcpp::stop("anchor_coords must have the same dimensions as coords"); anchor = matrix_to_points(a); }
    std::vector<double> anchorSchedule = resolve_schedule(anchor_weights, max_iter), bendSchedule = resolve_schedule(bend_weights, max_iter);
    std::vector<int> bendA = Rcpp::as<std::vector<int>>(bend_a), bendB = Rcpp::as<std::vector<int>>(bend_b), bendC = Rcpp::as<std::vector<int>>(bend_c);
    std::vector<int> trace_iteration; std::vector<double> trace_energy, trace_gmds_energy, trace_anchor_energy, trace_bend_energy, trace_gradient_norm, trace_step, trace_anchor_weight, trace_bend_weight; std::vector<bool> trace_accepted; Rcpp::List frames(return_trace ? max_iter + 1 : 1);
    const double eps2 = edge_length_epsilon * edge_length_epsilon, gradTol2 = grad_tol * grad_tol;

    BendState state = evaluate_state(current, cache, eps2, useAnchor ? &anchor : nullptr, anchorSchedule[0], bendA, bendB, bendC, bendSchedule[0]);
    trace_iteration.push_back(0); trace_energy.push_back(state.energy); trace_gmds_energy.push_back(state.gmdsEnergy); trace_anchor_energy.push_back(state.anchorEnergy); trace_bend_energy.push_back(state.bendEnergy); trace_gradient_norm.push_back(std::sqrt(state.gradNorm2)); trace_step.push_back(NA_REAL); trace_accepted.push_back(true); trace_anchor_weight.push_back(anchorSchedule[0]); trace_bend_weight.push_back(bendSchedule[0]); frames[0] = points_to_matrix(current, coords.ncol());

    int frameCount = 1;
    for(int iter = 1; iter <= max_iter; iter++){
        double iterAnchor = anchorSchedule[static_cast<size_t>(iter)], iterBend = bendSchedule[static_cast<size_t>(iter)];
        state = evaluate_state(current, cache, eps2, useAnchor ? &anchor : nullptr, iterAnchor, bendA, bendB, bendC, iterBend);
        if(!std::isfinite(state.energy) || state.gradNorm2 <= gradTol2) break;
        double step = initial_step; bool accepted = false; std::vector<Point<>> proposal = current; BendState candidate = state;
        while(std::isfinite(step) && step >= min_step){
            proposal = current;
            for(size_t i = 0; i < proposal.size(); i++) proposal[i] -= state.gradient[i] * step;
            if(recenter) recenter_points(proposal);
            candidate = evaluate_state(proposal, cache, eps2, useAnchor ? &anchor : nullptr, iterAnchor, bendA, bendB, bendC, iterBend);
            double target = state.energy - armijo_factor * step * state.gradNorm2;
            if(std::isfinite(candidate.energy) && candidate.energy <= target){ accepted = true; break; }
            step *= step_shrink;
        }
        trace_iteration.push_back(iter); trace_energy.push_back(accepted ? candidate.energy : state.energy); trace_gmds_energy.push_back(accepted ? candidate.gmdsEnergy : state.gmdsEnergy); trace_anchor_energy.push_back(accepted ? candidate.anchorEnergy : state.anchorEnergy); trace_bend_energy.push_back(accepted ? candidate.bendEnergy : state.bendEnergy); trace_gradient_norm.push_back(std::sqrt(accepted ? candidate.gradNorm2 : state.gradNorm2)); trace_step.push_back(accepted ? step : NA_REAL); trace_accepted.push_back(accepted); trace_anchor_weight.push_back(iterAnchor); trace_bend_weight.push_back(iterBend);
        if(!accepted) break;
        current = proposal;
        if(return_trace) frames[frameCount++] = points_to_matrix(current, coords.ncol());
    }
    if(!return_trace) {
        frames[0] = points_to_matrix(current, coords.ncol());
    } else if(frameCount < frames.size()) {
        Rcpp::List trimmed(frameCount);
        for(int i = 0; i < frameCount; i++) trimmed[i] = frames[i];
        frames = trimmed;
    }
    return Rcpp::List::create(Rcpp::_["coords"] = points_to_matrix(current, coords.ncol()),
                              Rcpp::_["trace"] = build_trace_df(trace_iteration, trace_energy, trace_gmds_energy, trace_anchor_energy, trace_bend_energy, trace_gradient_norm, trace_step, trace_accepted, trace_anchor_weight, trace_bend_weight),
                              Rcpp::_["frames"] = frames,
                              Rcpp::_["final_anchor_weight"] = trace_anchor_weight.back(),
                              Rcpp::_["final_bending_weight"] = trace_bend_weight.back(),
                              Rcpp::_["n_threads_used"] = 1L);
}
