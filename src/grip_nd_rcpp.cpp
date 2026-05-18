#include <Rcpp.h>

#include <cmath>
#include <ctime>
#include <limits>

#include "DrawGraphND.h"
#include "GraphND.h"

namespace {

void validate_weighted_nd_tuning(int dim,
                                 int rounds,
                                 int final_rounds,
                                 int num_init,
                                 int num_nbrs,
                                 double r,
                                 double s,
                                 double repulsion_factor,
                                 int tinit_factor)
{
    if(dim < 2)
        Rcpp::stop("dim must be >= 2");
    if(rounds < 0)
        Rcpp::stop("rounds must be a non-negative integer");
    if(final_rounds < 0)
        Rcpp::stop("final_rounds must be a non-negative integer");
    if(num_init < dim + 1)
        Rcpp::stop("num_init must be at least dim + 1 for ND weighted layout");
    if(num_nbrs <= 0)
        Rcpp::stop("num_nbrs must be a positive integer");
    if(!std::isfinite(r) || r < 0.0 || r > 1.0)
        Rcpp::stop("r must be finite and in [0, 1]");
    if(!std::isfinite(s) || s < 0.0)
        Rcpp::stop("s must be finite and >= 0");
    if(!std::isfinite(repulsion_factor) || repulsion_factor < 0.0)
        Rcpp::stop("repulsion_factor must be finite and >= 0");
    if(tinit_factor <= 0)
        Rcpp::stop("tinit_factor must be a positive integer");
}

gripnd::GraphND graph_from_r_lists(Rcpp::List adj_list,
                                   Rcpp::List weight_list,
                                   int n)
{
    if(adj_list.size() != n)
        Rcpp::stop("adj_list length must equal n");
    if(weight_list.size() != n)
        Rcpp::stop("weight_list length must equal n");

    gripnd::GraphND graph(static_cast<std::size_t>(n));
    for(int i = 0; i < n; i++){
        Rcpp::IntegerVector nbr = adj_list[i];
        Rcpp::NumericVector wt = weight_list[i];
        if(nbr.size() != wt.size())
            Rcpp::stop("weight_list entries must be parallel to adj_list entries");

        for(int k = 0; k < nbr.size(); k++){
            const int v = nbr[k];
            const double w = wt[k];
            if(v < 1 || v > n)
                Rcpp::stop("adj_list contains a vertex outside [1, n]");
            if(v == i + 1)
                continue;
            if(!std::isfinite(w) || w <= 0.0)
                Rcpp::stop("weight_list must contain finite values > 0");
            graph.add_directed_edge(
                static_cast<gripnd::vertex_t>(i),
                static_cast<gripnd::vertex_t>(v - 1),
                w
            );
        }
    }
    return graph;
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix grip_layout_weighted_nd_adj_cpp(Rcpp::List adj_list,
                                                    Rcpp::List weight_list,
                                                    int n,
                                                    int dim,
                                                    int rounds,
                                                    int final_rounds,
                                                    int num_init,
                                                    int num_nbrs,
                                                    double r,
                                                    double s,
                                                    double repulsion_factor,
                                                    int tinit_factor,
                                                    Rcpp::Nullable<int> seed = R_NilValue)
{
    if(n <= 0)
        Rcpp::stop("n must be a positive integer");
    validate_weighted_nd_tuning(
        dim,
        rounds,
        final_rounds,
        num_init,
        num_nbrs,
        r,
        s,
        repulsion_factor,
        tinit_factor
    );

    const unsigned int seed_value = seed.isNotNull()
        ? static_cast<unsigned int>(Rcpp::as<int>(seed))
        : static_cast<unsigned int>(std::time(NULL));

    gripnd::GraphND graph = graph_from_r_lists(adj_list, weight_list, n);
    gripnd::DrawGraphND drawer(
        graph,
        dim,
        rounds,
        final_rounds,
        num_init,
        num_nbrs,
        r,
        s,
        repulsion_factor,
        tinit_factor,
        seed_value
    );

    std::vector<gripnd::PointND> coords = drawer.layout();
    Rcpp::NumericMatrix out(n, dim);
    for(int i = 0; i < n; i++){
        for(int d = 0; d < dim; d++)
            out(i, d) = coords[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
    }
    return out;
}
