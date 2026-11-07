#include <Rcpp.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>

namespace {
using V = std::array<double, 3>;
V sub(const V& a, const V& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
double dot(const V& a, const V& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
V cross(const V& a, const V& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
double segment(const V& p, const V& a, const V& b) {
    V v=sub(b,a), w=sub(p,a);
    double vv=dot(v,v), t=vv > 0 ? std::max(0.0,std::min(1.0,dot(w,v)/vv)) : 0;
    for(int j=0;j<3;++j) w[j]-=t*v[j];
    return dot(w,w);
}
struct Face { V a,b,c,lo,hi,center; };
double distance2(const V& p, const Face& f) {
    V ab=sub(f.b,f.a), ac=sub(f.c,f.a), n=cross(ab,ac), ap=sub(p,f.a);
    double nn=dot(n,n);
    if(nn>0) {
        double t=dot(ap,n)/nn;
        V q=p;
        for(int j=0;j<3;++j) q[j]-=t*n[j];
        if(dot(cross(ab,sub(q,f.a)),n)>=0 &&
           dot(cross(sub(f.c,f.b),sub(q,f.b)),n)>=0 &&
           dot(cross(sub(f.a,f.c),sub(q,f.c)),n)>=0)
            return t*t*nn;
    }
    return std::min({segment(p,f.a,f.b),segment(p,f.b,f.c),segment(p,f.c,f.a)});
}
struct Node { V lo,hi; int begin,end,left=-1,right=-1; };
class Tree {
public:
    std::vector<Face> faces;
    std::vector<int> order;
    std::vector<Node> nodes;
    int build(int begin,int end) {
        Node node; node.begin=begin; node.end=end;
        node.lo=faces[order[begin]].lo; node.hi=faces[order[begin]].hi;
        for(int i=begin+1;i<end;++i) for(int j=0;j<3;++j) {
            node.lo[j]=std::min(node.lo[j],faces[order[i]].lo[j]);
            node.hi[j]=std::max(node.hi[j],faces[order[i]].hi[j]);
        }
        int index=nodes.size(); nodes.push_back(node);
        if(end-begin>8) {
            int axis=0;
            for(int j=1;j<3;++j) if(node.hi[j]-node.lo[j]>node.hi[axis]-node.lo[axis]) axis=j;
            int mid=begin+(end-begin)/2;
            std::nth_element(order.begin()+begin,order.begin()+mid,order.begin()+end,
                [&](int a,int b){return faces[a].center[axis]<faces[b].center[axis];});
            int l=build(begin,mid), r=build(mid,end);
            nodes[index].left=l; nodes[index].right=r;
        }
        return index;
    }
    double bound(const V& p,int id) const {
        double out=0;
        for(int j=0;j<3;++j) {
            double d=std::max({nodes[id].lo[j]-p[j],0.0,p[j]-nodes[id].hi[j]});
            out+=d*d;
        }
        return out;
    }
    void query(const V& p,int id,double& best) const {
        if(bound(p,id)>best) return;
        const Node& node=nodes[id];
        if(node.left<0) {
            for(int i=node.begin;i<node.end;++i) best=std::min(best,distance2(p,faces[order[i]]));
        } else {
            int a=node.left,b=node.right;
            if(bound(p,b)<bound(p,a)) std::swap(a,b);
            query(p,a,best); query(p,b,best);
        }
    }
};
}

// Exact closest-point distances to the union of triangles; the tree only prunes.
// [[Rcpp::export]]
Rcpp::NumericVector grip_surface_distances_cpp(Rcpp::NumericMatrix points,
        Rcpp::NumericMatrix vertices,Rcpp::IntegerMatrix triangles) {
    if(points.ncol()!=3 || vertices.ncol()!=3 || triangles.ncol()!=3 || triangles.nrow()<1)
        Rcpp::stop("Invalid surface distance dimensions");
    Tree tree;
    for(int i=0;i<triangles.nrow();++i) {
        if(i%4096==0) Rcpp::checkUserInterrupt();
        Face f;
        V* v[3]={&f.a,&f.b,&f.c};
        for(int k=0;k<3;++k) {
            int row=triangles(i,k)-1;
            if(row<0 || row>=vertices.nrow()) Rcpp::stop("Invalid triangle index");
            for(int j=0;j<3;++j) (*v[k])[j]=vertices(row,j);
        }
        for(int j=0;j<3;++j) {
            f.lo[j]=std::min({f.a[j],f.b[j],f.c[j]});
            f.hi[j]=std::max({f.a[j],f.b[j],f.c[j]});
            f.center[j]=(f.a[j]+f.b[j]+f.c[j])/3;
        }
        tree.faces.push_back(f);
    }
    tree.order.resize(tree.faces.size());
    std::iota(tree.order.begin(),tree.order.end(),0);
    tree.build(0,tree.faces.size());
    Rcpp::NumericVector out(points.nrow());
    for(int i=0;i<points.nrow();++i) {
        if(i%1024==0) Rcpp::checkUserInterrupt();
        V p={points(i,0),points(i,1),points(i,2)};
        double best=std::numeric_limits<double>::infinity();
        tree.query(p,0,best); out[i]=std::sqrt(best);
    }
    return out;
}
