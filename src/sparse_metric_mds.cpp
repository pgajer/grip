// Sparse constraints from Ortmann et al. (2017) / Zheng et al. (2018).
// Algorithm and pinned reference: dev/design/sparse-metric-mds.md.
// Adapted algorithm; MIT attribution retained in inst/COPYRIGHTS.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <vector>

namespace {
struct SparseTerm { int i, j; double d, a, b; };
using Link = std::pair<int, double>;
using Entry = std::pair<double, int>;
std::vector<double> distances(const std::vector<std::vector<Link>>& graph, int source) {
    const int n = graph.size();
    std::vector<double> d(n, std::numeric_limits<double>::infinity());
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    d[source] = 0; queue.push({0, source});
    std::size_t work = 0;
    while (!queue.empty()) {
        if ((work++ & 4095) == 0) Rcpp::checkUserInterrupt();
        const Entry e = queue.top(); queue.pop();
        if (e.first != d[e.second]) continue;
        for (const Link& edge : graph[e.second]) {
            const double candidate = e.first + edge.second;
            if (candidate < d[edge.first]) {
                d[edge.first] = candidate;
                queue.push({candidate, edge.first});
            }
        }
    }
    for (double value : d) if (!std::isfinite(value))
        Rcpp::stop("Sparse SGD requires a connected graph with representable path lengths");
    return d;
}
std::uint64_t sample_index(std::mt19937_64& rng, std::uint64_t bound) {
    const std::uint64_t threshold = -bound % bound;
    std::uint64_t x; do { x = rng(); } while (x < threshold);
    return x % bound;
}
}

// Returns original-unit targets and endpoint multiplicities, not dense APSP.
// [[Rcpp::export]]
Rcpp::List grip_sparse_prepare_cpp(int n, Rcpp::IntegerMatrix edges,
    Rcpp::NumericVector lengths, int h, Rcpp::IntegerVector supplied, int seed,
    double max_workspace_bytes, Rcpp::Nullable<Rcpp::NumericMatrix> observations = R_NilValue,
    std::string selection = "randomized", bool region_weighting = true,
    bool save_distances = false, std::string metric = "euclidean") {
    if (n < 2 || h < 1 || h > n || seed < 0 || edges.ncol() != 2 ||
        lengths.size() != edges.nrow() || (supplied.size() && supplied.size() != h) ||
        !std::isfinite(max_workspace_bytes) || max_workspace_bytes <= 0)
        Rcpp::stop("Invalid sparse preparation controls");
    const bool observation_input = observations.isNotNull();
    if (metric != "euclidean" && metric != "linf") Rcpp::stop("Unknown observation metric");
    if (metric == "linf" && !observation_input) Rcpp::stop("linf requires observations");
    Rcpp::NumericMatrix X;
    if (observation_input) {
        X = Rcpp::NumericMatrix(observations);
        if (X.nrow()!=n || X.ncol()<1) Rcpp::stop("Invalid observation dimensions");
        for (double z : X) if (!std::isfinite(z)) Rcpp::stop("Nonfinite observation");
    }
    if (selection!="randomized" && selection!="farthest" && selection!="uniform")
        Rcpp::stop("Unknown pivot selection rule");
    const std::size_t m = edges.nrow();
    const long double nh = static_cast<long double>(n) * h;
    // Conservative bound: distances, graph/queue, reserved terms, output, scratch.
    const long double estimate = (save_distances ? 16 : 8) * nh + 128 * m + 96 * n +
        96 * (nh + m);
    if (estimate > max_workspace_bytes || nh + m > std::numeric_limits<int>::max())
        Rcpp::stop("Sparse preparation exceeds max_workspace_bytes or constraint index range");
    if (m == 0 && !observation_input) Rcpp::stop("Sparse SGD requires a connected graph");
    double unit = m ? *std::max_element(lengths.begin(),lengths.end()) : 0;
    if (observation_input) for (double z : X) unit=std::max(unit,std::abs(z));
    if (!std::isfinite(unit) || unit<=0) Rcpp::stop("Invalid distance scale");
    std::vector<std::vector<Link>> graph(n);
    for (std::size_t e = 0; e < m; ++e) {
        if (edges(e,0) == NA_INTEGER || edges(e,1) == NA_INTEGER) Rcpp::stop("Invalid sparse graph edge");
        int i = edges(e,0)-1, j = edges(e,1)-1;
        if (i < 0 || j < 0 || i >= n || j >= n || i == j ||
            !std::isfinite(lengths[e]) || lengths[e] <= 0)
            Rcpp::stop("Invalid sparse graph edge");
        const double length = lengths[e] / unit;
        if (!std::isfinite(length) || length <= 0) Rcpp::stop("Normalized edge lengths exceed numeric range");
        graph[i].push_back({j,length}); graph[j].push_back({i,length});
    }
    for (auto& row : graph) {
        std::sort(row.begin(), row.end());
        for (std::size_t j=1; j<row.size(); ++j)
            if (row[j-1].first == row[j].first) Rcpp::stop("Sparse graph edges must be unique");
    }
    std::vector<bool> selected(n, false);
    for (int p : supplied) {
        if (p < 1 || p > n || selected[p-1]) Rcpp::stop("pivots must be distinct valid vertex ids");
        selected[p-1] = true;
    }
    std::fill(selected.begin(), selected.end(), false);
    std::vector<int> pivots(h), owner(n, -1);
    std::vector<double> nearest(n, std::numeric_limits<double>::infinity());
    std::vector<std::vector<double>> all; all.reserve(h);
    std::mt19937_64 rng(seed);
    for (int k=0; k<h; ++k) {
        Rcpp::checkUserInterrupt();
        int p = -1;
        if (supplied.size()) p = supplied[k]-1;
        else if (k == 0) p = sample_index(rng, n);
        else if (selection == "uniform") {
            int rank = sample_index(rng,n-k);
            for (int i=0; i<n; ++i) if (!selected[i] && rank-- == 0) { p=i; break; }
        } else if (selection == "farthest") {
            for (int i=0; i<n; ++i) if (!selected[i] && (p<0 || nearest[i]>nearest[p])) p=i;
        } else {
            double largest = *std::max_element(nearest.begin(), nearest.end());
            long double total = 0;
            for (int i=0; i<n; ++i) if (!selected[i]) total += nearest[i]/largest;
            long double draw = static_cast<double>(rng() >> 11) * 0x1.0p-53 * total;
            long double cumulative = 0;
            for (int i=0; i<n; ++i) if (!selected[i]) {
                p = i; cumulative += nearest[i]/largest;
                if (draw < cumulative) break;
            }
        }
        if (p < 0 || selected[p]) Rcpp::stop("Sparse pivot sampling failed");
        selected[p] = true; pivots[k] = p;
        if (observation_input) {
            std::vector<double> row(n,0);
            if (metric == "linf") {
                // R matrices are column-major: stream each feature contiguously.
                for (int col=0; col<X.ncol(); ++col) {
                    const double pivot_value = X(p,col);
                    for (int i=0; i<n; ++i)
                        row[i] = std::max(row[i],std::abs((X(i,col)-pivot_value)/unit));
                }
            } else {
                for (int i=0; i<n; ++i) if (i!=p)
                    for (int col=0; col<X.ncol(); ++col)
                        row[i]=std::hypot(row[i],(X(i,col)-X(p,col))/unit);
            }
            for (int i=0; i<n; ++i) if (i!=p && (!std::isfinite(row[i]) || row[i]<=0))
                Rcpp::stop("Observation distances must be positive and representable; consolidate duplicates");
            all.push_back(std::move(row));
        } else all.push_back(distances(graph,p));
        for (int i=0; i<n; ++i) if (all.back()[i] < nearest[i]) {
            nearest[i] = all.back()[i]; owner[i] = k;
        }
    }
    std::vector<std::vector<double>> region(h);
    for (int i=0; i<n; ++i) region[owner[i]].push_back(nearest[i]);
    for (auto& row : region) std::sort(row.begin(),row.end());
    std::vector<SparseTerm> terms;
    terms.reserve(static_cast<std::size_t>(nh)+m);
    for (std::size_t e=0; e<m; ++e) {
        int i = edges(e,0)-1, j = edges(e,1)-1;
        if (i>j) std::swap(i,j);
        terms.push_back({i,j,lengths[e],1,1});
    }
    for (int k=0; k<h; ++k) {
        Rcpp::checkUserInterrupt();
        const int p=pivots[k];
        for (int i=0; i<n; ++i) {
            if (i == p) continue;
            const auto found = std::lower_bound(graph[p].begin(),graph[p].end(),Link(i,0));
            if (found != graph[p].end() && found->first == i) continue;
            const double d = all[k][i];
            const double s = region_weighting ?
                std::upper_bound(region[k].begin(),region[k].end(),d/2)-region[k].begin() : 1;
            if (!region_weighting && selected[i] && i<p) continue;
            const double target = d * unit;
            if (!std::isfinite(target) || target <= 0) Rcpp::stop("Sparse target distances exceed numeric range");
            if (!region_weighting) terms.push_back({std::min(i,p),std::max(i,p),target,1,1});
            else if (i<p) terms.push_back({i,p,target,s,0});
            else terms.push_back({p,i,target,0,s});
        }
    }
    std::sort(terms.begin(), terms.end(), [](const SparseTerm& a,const SparseTerm& b) {
        return a.i < b.i || (a.i == b.i && a.j < b.j);
    });
    std::size_t count=0;
    for (std::size_t k=0; k<terms.size(); ++k) {
        if ((k & 16383) == 0) Rcpp::checkUserInterrupt();
        if (count && terms[count-1].i==terms[k].i && terms[count-1].j==terms[k].j) {
            // Opposite Dijkstra directions can round differently; pick min.
            terms[count-1].d = std::min(terms[count-1].d, terms[k].d);
            terms[count-1].a += terms[k].a; terms[count-1].b += terms[k].b;
        } else terms[count++] = terms[k];
    }
    Rcpp::IntegerMatrix pairs(count,2);
    Rcpp::NumericVector target(count), ci(count), cj(count);
    for (std::size_t k=0; k<count; ++k) {
        pairs(k,0)=terms[k].i+1; pairs(k,1)=terms[k].j+1;
        target[k]=terms[k].d; ci[k]=terms[k].a; cj[k]=terms[k].b;
    }
    Rcpp::IntegerVector p(h), r(n);
    for (int k=0; k<h; ++k) p[k]=pivots[k]+1;
    for (int i=0; i<n; ++i) r[i]=pivots[owner[i]]+1;
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("pairs")=pairs,Rcpp::Named("targets")=target,
        Rcpp::Named("count_i")=ci,Rcpp::Named("count_j")=cj,
        Rcpp::Named("pivots")=p,Rcpp::Named("region")=r,
        Rcpp::Named("workspace_estimate_bytes")=static_cast<double>(estimate));
    if (save_distances) {
        Rcpp::NumericMatrix rows(h,n);
        for (int k=0; k<h; ++k) for (int i=0; i<n; ++i) rows(k,i)=all[k][i]*unit;
        out["distance_rows"]=rows;
    }
    return out;
}

// Exact L-infinity kNN, with O(n*p + n*q) storage and deterministic ties.
// [[Rcpp::export]]
Rcpp::List grip_linf_knn_cpp(Rcpp::NumericMatrix X, int q, double max_workspace_bytes) {
    const int n = X.nrow(), p = X.ncol();
    if (n < 2 || p < 1 || q < 1 || q >= n ||
        !std::isfinite(max_workspace_bytes) || max_workspace_bytes <= 0)
        Rcpp::stop("Invalid L-infinity neighbor controls");
    const long double nq = static_cast<long double>(n)*q;
    const long double estimate = 8*static_cast<long double>(n)*p + 64*nq + 32*n;
    if (estimate > max_workspace_bytes || nq > std::numeric_limits<int>::max())
        Rcpp::stop("L-infinity neighbor search exceeds workspace or edge index range");
    // Contiguous observations allow early rejection once a coordinate exceeds
    // the current kth distance. This pruning never changes the exact result.
    std::vector<double> rows(static_cast<std::size_t>(n)*p);
    for (int col=0; col<p; ++col) for (int i=0; i<n; ++i) {
        const double value = X(i,col);
        if (!std::isfinite(value)) Rcpp::stop("Nonfinite observation");
        rows[static_cast<std::size_t>(i)*p+col] = value;
    }
    std::vector<SparseTerm> edges; edges.reserve(static_cast<std::size_t>(nq));
    for (int i=0; i<n; ++i) {
        Rcpp::checkUserInterrupt();
        std::priority_queue<Entry> nearest;
        const double* a = &rows[static_cast<std::size_t>(i)*p];
        for (int j=0; j<n; ++j) if (i != j) {
            if ((j & 4095) == 0) Rcpp::checkUserInterrupt();
            const double* b = &rows[static_cast<std::size_t>(j)*p];
            double d = 0;
            const double limit = nearest.size() == static_cast<std::size_t>(q) ?
                nearest.top().first : std::numeric_limits<double>::infinity();
            for (int col=0; col<p; ++col) {
                d = std::max(d,std::abs(a[col]-b[col]));
                if (d > limit) break;
            }
            if (d == 0) Rcpp::stop("Duplicate observations in L-infinity neighbor search");
            const Entry candidate(d,j);
            if (nearest.size() < static_cast<std::size_t>(q)) nearest.push(candidate);
            else if (candidate < nearest.top()) { nearest.pop(); nearest.push(candidate); }
        }
        while (!nearest.empty()) {
            const Entry e = nearest.top(); nearest.pop();
            if (!std::isfinite(e.first)) Rcpp::stop("L-infinity edge distances exceed numeric range");
            edges.push_back({std::min(i,e.second),std::max(i,e.second),e.first,1,1});
        }
    }
    std::sort(edges.begin(),edges.end(),[](const SparseTerm& a,const SparseTerm& b) {
        return a.i < b.i || (a.i == b.i && a.j < b.j);
    });
    edges.erase(std::unique(edges.begin(),edges.end(),[](const SparseTerm& a,const SparseTerm& b) {
        return a.i == b.i && a.j == b.j;
    }),edges.end());
    Rcpp::IntegerMatrix pairs(edges.size(),2);
    Rcpp::NumericVector lengths(edges.size());
    for (std::size_t k=0; k<edges.size(); ++k) {
        pairs(k,0)=edges[k].i+1; pairs(k,1)=edges[k].j+1; lengths[k]=edges[k].d;
    }
    return Rcpp::List::create(Rcpp::Named("pairs")=pairs,Rcpp::Named("lengths")=lengths);
}
