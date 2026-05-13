// C++ backend for edge-isometric repulsive unfolding.
// The R implementation remains the behavioral reference.

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

struct EdgeRepulsiveState {
  double energy;
  double edge_energy;
  double repel_energy;
  NumericMatrix gradient;
  double gradient_norm;
  bool feasible;
  NumericVector edge_lengths;
  NumericVector edge_relative_lengths;
  NumericVector edge_residuals;
  int n_wall_violations;
};

struct RepulsiveState {
  double energy;
  double repel_energy;
  NumericMatrix gradient;
  double gradient_norm;
  NumericVector pair_lengths;
};

static double sqr(double x) {
  return x * x;
}

static void center_coords(NumericMatrix& Z) {
  int n = Z.nrow();
  int d = Z.ncol();
  for (int j = 0; j < d; ++j) {
    double mean = 0.0;
    for (int i = 0; i < n; ++i) mean += Z(i, j);
    mean /= static_cast<double>(n);
    for (int i = 0; i < n; ++i) Z(i, j) -= mean;
  }
}

static EdgeRepulsiveState compute_state_cpp(const NumericMatrix& Z,
                                            const IntegerMatrix& edges,
                                            const NumericVector& edge_lengths,
                                            const NumericVector& edge_weights,
                                            const std::string& edge_family,
                                            double eps_plus,
                                            double beta,
                                            double lambda,
                                            const IntegerMatrix& pair_index,
                                            const NumericVector& pair_weights,
                                            const std::string& repulsion_family,
                                            double repulsion_delta,
                                            double repulsion_power,
                                            double distance_eps) {
  int n = Z.nrow();
  int dim = Z.ncol();
  int m = edges.nrow();
  EdgeRepulsiveState out;
  out.energy = 0.0;
  out.edge_energy = 0.0;
  out.repel_energy = 0.0;
  out.gradient = NumericMatrix(n, dim);
  out.gradient_norm = 0.0;
  out.feasible = true;
  out.edge_lengths = NumericVector(m);
  out.edge_relative_lengths = NumericVector(m);
  out.edge_residuals = NumericVector(m);
  out.n_wall_violations = 0;

  for (int e = 0; e < m; ++e) {
    int u = edges(e, 0) - 1;
    int v = edges(e, 1) - 1;
    double dist2 = distance_eps * distance_eps;
    for (int j = 0; j < dim; ++j) {
      double diff = Z(u, j) - Z(v, j);
      dist2 += diff * diff;
    }
    double dist = std::sqrt(dist2);
    double ell = edge_lengths[e];
    double s = dist / ell;
    double value = 0.0;
    double deriv = 0.0;

    if (edge_family == "upper_barrier" && beta > 0.0) {
      double uu = (s - 1.0) / eps_plus;
      if (!(uu < 1.0 - 1e-10) || !R_finite(uu)) {
        out.feasible = false;
        out.n_wall_violations += 1;
        value = R_PosInf;
        deriv = NA_REAL;
      } else {
        value = 0.5 * sqr(s - 1.0) + beta * (-std::log1p(-uu) - uu);
        deriv = (s - 1.0) + beta * (1.0 / (1.0 - uu) - 1.0) / eps_plus;
      }
    } else {
      value = 0.5 * sqr(s - 1.0);
      deriv = s - 1.0;
    }

    out.edge_lengths[e] = dist;
    out.edge_relative_lengths[e] = s;
    out.edge_residuals[e] = dist - ell;
    if (R_finite(value)) {
      out.edge_energy += edge_weights[e] * value;
      double coeff = edge_weights[e] * deriv / ell / dist;
      for (int j = 0; j < dim; ++j) {
        double g = coeff * (Z(u, j) - Z(v, j));
        out.gradient(u, j) += g;
        out.gradient(v, j) -= g;
      }
    } else {
      out.edge_energy = R_PosInf;
    }
  }

  // Match the R reference: if any edge barrier is infeasible, the edge term
  // contributes infinite energy and no partial edge gradient is reported.
  if (!R_finite(out.edge_energy)) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < dim; ++j) {
        out.gradient(i, j) = 0.0;
      }
    }
  }

  if (lambda > 0.0) {
    int np = pair_index.nrow();
    for (int p = 0; p < np; ++p) {
      int u = pair_index(p, 0) - 1;
      int v = pair_index(p, 1) - 1;
      double dist2 = distance_eps * distance_eps;
      for (int j = 0; j < dim; ++j) {
        double diff = Z(u, j) - Z(v, j);
        dist2 += diff * diff;
      }
      double dist = std::sqrt(dist2);
      double rr = dist * dist + repulsion_delta * repulsion_delta;
      double value = 0.0;
      double deriv_r = 0.0;
      if (repulsion_family == "inverse_power") {
        value = std::pow(rr, -repulsion_power / 2.0);
        deriv_r = -repulsion_power * dist * std::pow(rr, -repulsion_power / 2.0 - 1.0);
      } else {
        value = -0.5 * std::log(rr);
        deriv_r = -dist / rr;
      }
      double w = pair_weights[p];
      out.repel_energy += w * value;
      double coeff = lambda * w * deriv_r / dist;
      for (int j = 0; j < dim; ++j) {
        double g = coeff * (Z(u, j) - Z(v, j));
        out.gradient(u, j) += g;
        out.gradient(v, j) -= g;
      }
    }
  }

  out.energy = out.edge_energy + lambda * out.repel_energy;
  double gn = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < dim; ++j) {
      gn += out.gradient(i, j) * out.gradient(i, j);
    }
  }
  out.gradient_norm = std::sqrt(gn);
  return out;
}

static RepulsiveState compute_repulsive_state_cpp(const NumericMatrix& Z,
                                                  double lambda,
                                                  const IntegerMatrix& pair_index,
                                                  const NumericVector& pair_weights,
                                                  const std::string& repulsion_family,
                                                  double repulsion_delta,
                                                  double repulsion_power,
                                                  double distance_eps) {
  int n = Z.nrow();
  int dim = Z.ncol();
  int np = pair_index.nrow();
  RepulsiveState out;
  out.energy = 0.0;
  out.repel_energy = 0.0;
  out.gradient = NumericMatrix(n, dim);
  out.gradient_norm = 0.0;
  out.pair_lengths = NumericVector(np);

  if (lambda > 0.0) {
    for (int p = 0; p < np; ++p) {
      int u = pair_index(p, 0) - 1;
      int v = pair_index(p, 1) - 1;
      double dist2 = distance_eps * distance_eps;
      for (int j = 0; j < dim; ++j) {
        double diff = Z(u, j) - Z(v, j);
        dist2 += diff * diff;
      }
      double dist = std::sqrt(dist2);
      out.pair_lengths[p] = dist;
      double rr = dist * dist + repulsion_delta * repulsion_delta;
      double value = 0.0;
      double deriv_r = 0.0;
      if (repulsion_family == "inverse_power") {
        value = std::pow(rr, -repulsion_power / 2.0);
        deriv_r = -repulsion_power * dist * std::pow(rr, -repulsion_power / 2.0 - 1.0);
      } else {
        value = -0.5 * std::log(rr);
        deriv_r = -dist / rr;
      }
      double w = pair_weights[p];
      out.repel_energy += w * value;
      double coeff = lambda * w * deriv_r / dist;
      for (int j = 0; j < dim; ++j) {
        double g = coeff * (Z(u, j) - Z(v, j));
        out.gradient(u, j) += g;
        out.gradient(v, j) -= g;
      }
    }
  } else {
    for (int p = 0; p < np; ++p) {
      int u = pair_index(p, 0) - 1;
      int v = pair_index(p, 1) - 1;
      double dist2 = distance_eps * distance_eps;
      for (int j = 0; j < dim; ++j) {
        double diff = Z(u, j) - Z(v, j);
        dist2 += diff * diff;
      }
      out.pair_lengths[p] = std::sqrt(dist2);
    }
  }

  out.energy = lambda * out.repel_energy;
  double gn = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < dim; ++j) {
      gn += out.gradient(i, j) * out.gradient(i, j);
    }
  }
  out.gradient_norm = std::sqrt(gn);
  return out;
}

static List state_to_list(const EdgeRepulsiveState& s) {
  return List::create(
    _["energy"] = s.energy,
    _["edge.energy"] = s.edge_energy,
    _["repel.energy"] = s.repel_energy,
    _["gradient"] = s.gradient,
    _["gradient.norm"] = s.gradient_norm,
    _["feasible"] = s.feasible,
    _["edge.embedded.lengths"] = s.edge_lengths,
    _["edge.relative.lengths"] = s.edge_relative_lengths,
    _["edge.residuals"] = s.edge_residuals,
    _["n.wall.violations"] = s.n_wall_violations
  );
}

static List repulsive_state_to_list(const RepulsiveState& s) {
  return List::create(
    _["energy"] = s.energy,
    _["repel.energy"] = s.repel_energy,
    _["gradient"] = s.gradient,
    _["gradient.norm"] = s.gradient_norm,
    _["pair.embedded.lengths"] = s.pair_lengths
  );
}

// [[Rcpp::export]]
List grip_edge_repulsive_state_cpp(NumericMatrix coords,
                                   IntegerMatrix edges,
                                   NumericVector edge_lengths,
                                   NumericVector edge_weights,
                                   std::string edge_family,
                                   double eps_plus,
                                   double beta,
                                   double lambda,
                                   IntegerMatrix pair_index,
                                   NumericVector pair_weights,
                                   std::string repulsion_family,
                                   double repulsion_delta,
                                   double repulsion_power,
                                   double distance_eps) {
  EdgeRepulsiveState s = compute_state_cpp(
    coords, edges, edge_lengths, edge_weights, edge_family, eps_plus, beta,
    lambda, pair_index, pair_weights, repulsion_family, repulsion_delta,
    repulsion_power, distance_eps
  );
  return state_to_list(s);
}

// [[Rcpp::export]]
List grip_repulsive_state_cpp(NumericMatrix coords,
                              double lambda,
                              IntegerMatrix pair_index,
                              NumericVector pair_weights,
                              std::string repulsion_family,
                              double repulsion_delta,
                              double repulsion_power,
                              double distance_eps) {
  RepulsiveState s = compute_repulsive_state_cpp(
    coords, lambda, pair_index, pair_weights, repulsion_family,
    repulsion_delta, repulsion_power, distance_eps
  );
  return repulsive_state_to_list(s);
}

// [[Rcpp::export]]
List grip_optimize_edge_repulsive_stage_cpp(NumericMatrix coords,
                                            IntegerMatrix edges,
                                            NumericVector edge_lengths,
                                            NumericVector edge_weights,
                                            double lambda,
                                            std::string edge_family,
                                            double eps_plus,
                                            double beta,
                                            IntegerMatrix pair_index,
                                            NumericVector pair_weights,
                                            std::string repulsion_family,
                                            double repulsion_delta,
                                            double repulsion_power,
                                            int max_iter,
                                            double initial_step,
                                            double step_shrink,
                                            double armijo,
                                            double min_step,
                                            double grad_tol,
                                            bool recenter,
                                            double distance_eps,
                                            bool return_frames) {
  NumericMatrix Z = clone(coords);
  if (recenter) center_coords(Z);
  EdgeRepulsiveState state = compute_state_cpp(
    Z, edges, edge_lengths, edge_weights, edge_family, eps_plus, beta,
    lambda, pair_index, pair_weights, repulsion_family, repulsion_delta,
    repulsion_power, distance_eps
  );

  std::vector<int> iteration;
  std::vector<double> energy, edge_energy, repel_energy, gradient_norm, step_vec;
  std::vector<int> accepted, wall, frame_index;
  List frames(return_frames ? max_iter + 1 : 0);
  int n_frames = 0;
  if (return_frames) {
    frames[n_frames++] = clone(Z);
  }

  auto push_trace = [&](int iter, const EdgeRepulsiveState& s, double step, bool accept, int frame) {
    iteration.push_back(iter);
    energy.push_back(s.energy);
    edge_energy.push_back(s.edge_energy);
    repel_energy.push_back(s.repel_energy);
    gradient_norm.push_back(s.gradient_norm);
    step_vec.push_back(step);
    accepted.push_back(accept ? 1 : 0);
    wall.push_back(s.n_wall_violations);
    frame_index.push_back(frame);
  };

  for (int iter = 0; iter < max_iter; ++iter) {
    push_trace(iter, state, NA_REAL, true, n_frames > 0 ? n_frames - 1 : NA_INTEGER);
    if (!R_finite(state.energy) || !R_finite(state.gradient_norm) ||
        state.gradient_norm <= grad_tol) {
      break;
    }
    double step = initial_step;
    bool did_accept = false;
    NumericMatrix candidate = clone(Z);
    EdgeRepulsiveState candidate_state = state;
    while (R_finite(step) && step >= min_step) {
      NumericMatrix proposal = clone(Z);
      for (int i = 0; i < proposal.nrow(); ++i) {
        for (int j = 0; j < proposal.ncol(); ++j) {
          proposal(i, j) -= step * state.gradient(i, j);
        }
      }
      if (recenter) center_coords(proposal);
      EdgeRepulsiveState proposal_state = compute_state_cpp(
        proposal, edges, edge_lengths, edge_weights, edge_family, eps_plus, beta,
        lambda, pair_index, pair_weights, repulsion_family, repulsion_delta,
        repulsion_power, distance_eps
      );
      double target = state.energy - armijo * step * state.gradient_norm * state.gradient_norm;
      if (proposal_state.feasible && R_finite(proposal_state.energy) &&
          proposal_state.energy <= target) {
        did_accept = true;
        candidate = proposal;
        candidate_state = proposal_state;
        break;
      }
      step *= step_shrink;
    }
    if (!did_accept) break;
    Z = candidate;
    state = candidate_state;
    if (return_frames) {
      frames[n_frames++] = clone(Z);
    }
    push_trace(iter + 1, state, step, true, n_frames > 0 ? n_frames - 1 : NA_INTEGER);
  }

  DataFrame trace = DataFrame::create(
    _["iteration"] = iteration,
    _["energy"] = energy,
    _["edge.energy"] = edge_energy,
    _["repel.energy"] = repel_energy,
    _["gradient.norm"] = gradient_norm,
    _["step"] = step_vec,
    _["accepted"] = accepted,
    _["n.wall.violations"] = wall,
    _["frame.index"] = frame_index
  );

  if (return_frames) {
    List trimmed(n_frames);
    for (int i = 0; i < n_frames; ++i) trimmed[i] = frames[i];
    frames = trimmed;
  }

  return List::create(
    _["coords"] = Z,
    _["state"] = state_to_list(state),
    _["trace"] = trace,
    _["frames"] = return_frames ? static_cast<SEXP>(frames) : R_NilValue
  );
}

// [[Rcpp::export]]
List grip_optimize_repulsive_stage_cpp(NumericMatrix coords,
                                       double lambda,
                                       IntegerMatrix pair_index,
                                       NumericVector pair_weights,
                                       std::string repulsion_family,
                                       double repulsion_delta,
                                       double repulsion_power,
                                       int max_iter,
                                       double initial_step,
                                       double step_shrink,
                                       double armijo,
                                       double min_step,
                                       double grad_tol,
                                       bool recenter,
                                       double distance_eps,
                                       bool return_frames) {
  NumericMatrix Z = clone(coords);
  if (recenter) center_coords(Z);
  RepulsiveState state = compute_repulsive_state_cpp(
    Z, lambda, pair_index, pair_weights, repulsion_family,
    repulsion_delta, repulsion_power, distance_eps
  );

  std::vector<int> iteration;
  std::vector<double> energy, repel_energy, gradient_norm, step_vec;
  std::vector<int> accepted, frame_index;
  List frames(return_frames ? max_iter + 1 : 0);
  int n_frames = 0;
  if (return_frames) {
    frames[n_frames++] = clone(Z);
  }

  auto push_trace = [&](int iter, const RepulsiveState& s, double step, bool accept, int frame) {
    iteration.push_back(iter);
    energy.push_back(s.energy);
    repel_energy.push_back(s.repel_energy);
    gradient_norm.push_back(s.gradient_norm);
    step_vec.push_back(step);
    accepted.push_back(accept ? 1 : 0);
    frame_index.push_back(frame);
  };

  for (int iter = 0; iter < max_iter; ++iter) {
    push_trace(iter, state, NA_REAL, true, n_frames > 0 ? n_frames - 1 : NA_INTEGER);
    if (!R_finite(state.energy) || !R_finite(state.gradient_norm) ||
        state.gradient_norm <= grad_tol) {
      break;
    }
    double step = initial_step;
    bool did_accept = false;
    NumericMatrix candidate = clone(Z);
    RepulsiveState candidate_state = state;
    while (R_finite(step) && step >= min_step) {
      NumericMatrix proposal = clone(Z);
      for (int i = 0; i < proposal.nrow(); ++i) {
        for (int j = 0; j < proposal.ncol(); ++j) {
          proposal(i, j) -= step * state.gradient(i, j);
        }
      }
      if (recenter) center_coords(proposal);
      RepulsiveState proposal_state = compute_repulsive_state_cpp(
        proposal, lambda, pair_index, pair_weights, repulsion_family,
        repulsion_delta, repulsion_power, distance_eps
      );
      double target = state.energy - armijo * step * state.gradient_norm * state.gradient_norm;
      if (R_finite(proposal_state.energy) && proposal_state.energy <= target) {
        did_accept = true;
        candidate = proposal;
        candidate_state = proposal_state;
        break;
      }
      step *= step_shrink;
    }
    if (!did_accept) break;
    Z = candidate;
    state = candidate_state;
    if (return_frames) {
      frames[n_frames++] = clone(Z);
    }
    push_trace(iter + 1, state, step, true, n_frames > 0 ? n_frames - 1 : NA_INTEGER);
  }

  DataFrame trace = DataFrame::create(
    _["iteration"] = iteration,
    _["energy"] = energy,
    _["repel.energy"] = repel_energy,
    _["gradient.norm"] = gradient_norm,
    _["step"] = step_vec,
    _["accepted"] = accepted,
    _["frame.index"] = frame_index
  );

  if (return_frames) {
    List trimmed(n_frames);
    for (int i = 0; i < n_frames; ++i) trimmed[i] = frames[i];
    frames = trimmed;
  }

  return List::create(
    _["coords"] = Z,
    _["state"] = repulsive_state_to_list(state),
    _["trace"] = trace,
    _["frames"] = return_frames ? static_cast<SEXP>(frames) : R_NilValue
  );
}
