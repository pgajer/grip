#!/usr/bin/env python3
"""Diagnostic checks of stated algebra and counterexamples, not proof certificates."""
import json
from pathlib import Path
import numpy as np

P = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(20260906)
checks = {}

# Recompute the centered tree-distance matrix independently of its factorization.
labels = np.repeat(np.arange(4), 3)
t = rng.uniform(.1, 1, len(labels))
F = np.eye(4)[labels] * t[:, None]
J = np.eye(len(t))-np.ones((len(t), len(t)))/len(t)
T = np.where(labels[:, None] == labels[None, :], abs(t[:, None]-t[None, :]), t[:, None]+t[None, :])
direct = -.5*J@(T*T)@J
factor = J@F@(2*np.eye(4)-np.ones((4,4)))@F.T@J
err = np.max(abs(direct-factor))
assert err < 1e-13
ev = np.linalg.eigvalsh(direct)
assert np.linalg.matrix_rank(J@F) == 4
assert sum(ev>1e-10) == 3 and sum(ev < -1e-10) == 1
checks['tree_factorization'] = dict(max_difference=float(err), positive=3, negative=1)

# Boundary sample designs: four equal tips, and the same tips plus the root.
tip = 2*(np.ones((4,4))-np.eye(4)); j4 = np.eye(4)-np.ones((4,4))/4
assert np.allclose(-.5*j4@(tip*tip)@j4, 2*j4)
root = np.ones((5,5)); root[1:,1:]=tip; root[0,0]=0
j5=np.eye(5)-np.ones((5,5))/5
evroot=np.linalg.eigvalsh(-.5*j5@(root*root)@j5)
assert sum(evroot>1e-10)==3 and sum(evroot < -1e-10)==1
checks['tree_boundary_designs']=dict(equal_tips_positive=3, root_and_tips_positive=3,
                                    root_and_tips_negative=1)

# Graph-metric acceleration gives the radial convexity used for disk containment.
errors=[];lower_gaps=[]
for eta in [-1,1]:
    for _ in range(100):
        xy=rng.normal(size=2);vel=rng.normal(size=2)
        gradient=2*xy*np.array([1,eta]);hessian=2*np.diag([1,eta])
        metric=np.eye(2)+np.outer(gradient,gradient)
        acceleration=-np.linalg.solve(metric,gradient)*(vel@hessian@vel)
        rho2=xy@xy;v2=vel@vel
        second=2*v2+2*xy@acceleration
        expression=2*v2-8*(xy[0]**2+eta*xy[1]**2)*(vel[0]**2+eta*vel[1]**2)/(1+4*rho2)
        errors.append(abs(second-expression));lower_gaps.append(second-2*v2/(1+4*rho2))
assert max(errors)<1e-12 and min(lower_gaps)>-1e-12
checks['disk_geodesic_containment']=dict(max_formula_difference=float(max(errors)),
                                       minimum_convexity_gap=float(min(lower_gaps)))

# Profile a target scale from weighted least squares at a generic nonzero shape.
d = rng.uniform(.2, 2, 20); ell = rng.uniform(.1, 1, 20); w = rng.uniform(.2, 3, 20)
def energy(d):
    a = np.sum(w*d*ell)/np.sum(w*ell**2)
    return .5*np.sum(w*(d-a*ell)**2), a
base, a = energy(d)
errs=[]
for c in [.01,.1,.5,2,10]:
    val, ac=energy(c*d)
    errs.extend([abs(val/(c*c)-base),abs(ac/c-a)])
A=np.sum(w*d*d);B=np.sum(w*ell*ell);C=np.sum(w*d*ell)
Q=1-C*C/(A*B);c=C/A
fixed=.5*np.sum(w*(c*d-ell)**2)
assert max(errs) < 1e-12 and abs(fixed-B*Q/2)<1e-12
checks['profile_homogeneity'] = dict(max_difference=float(max(errs)), fixed_scale_identity_difference=float(abs(fixed-B*Q/2)))

# Euler identity from an independently assembled coordinate gradient on a complete graph.
zz = rng.normal(size=(5,3)); ii,jj = np.triu_indices(5,1)
vv = zz[ii]-zz[jj]; dd = np.linalg.norm(vv,axis=1)
ll = rng.uniform(.2,2,len(dd)); ww = rng.uniform(.3,2,len(dd))
aa = np.sum(ww*dd*ll)/np.sum(ww*ll*ll)
ee = .5*np.sum(ww*(dd-aa*ll)**2)
gradient = np.zeros_like(zz)
edge_gradient = (ww*(dd-aa*ll)/dd)[:,None]*vv
np.add.at(gradient,ii,edge_gradient); np.add.at(gradient,jj,-edge_gradient)
euler_gap = abs(np.sum(zz*gradient)-2*ee)
assert ee>0 and euler_gap<1e-12
assert .5*np.sum(ww*(.99*dd-.99*aa*ll)**2) < ee
assert np.sum(ww*(np.zeros_like(dd)-0*ll)**2) == 0
checks['profile_euler_and_contraction'] = dict(euler_gap=float(euler_gap), positive_loss=float(ee))

# The removed historical update raises stress on an exact unit path realization.
z=np.array([-1.,0.,1.]); Wedge=np.array([[0.,2.,0.],[2.,0.,2.],[0.,2.,0.]])
V=np.diag(Wedge.sum(1))-Wedge; W=np.diag(abs(V).sum(1))
new=np.linalg.pinv(W)@V@z
lengths=np.array([abs(new[1]-new[0]),abs(new[2]-new[1]),abs(new[1]-new[0])+abs(new[2]-new[1])])
stress=np.sum((lengths-np.array([1,1,2]))**2)
assert np.allclose(new,[-.5,0,.5]) and stress==1.5
checks['removed_update_counterexample'] = dict(before_stress=0,after_stress=float(stress))

# Fixed chosen shortest paths can violate a triangle after coordinates move.
# Input: unit square cycle 0--1--2--3--0; choose 0--1--2 for the tied 0-to-2 route.
z=np.array([[0.,0.],[0.,10.],[2.,0.],[1.,0.]])
l02=np.linalg.norm(z[0]-z[1])+np.linalg.norm(z[1]-z[2])
l03=np.linalg.norm(z[0]-z[3]);l32=np.linalg.norm(z[3]-z[2])
assert l02 > l03+l32
checks['fixed_route_triangle_counterexample'] = dict(long_route=float(l02),two_direct_edges=float(l03+l32))

result=dict(status='passed',scope='Algebraic diagnostics and explicit counterexamples; analytic proofs require separate reading.',checks=checks)
(P/'build/validation').mkdir(parents=True, exist_ok=True)
(P/'build/validation/mathematics-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(f'All {len(checks)} algebra/counterexample checks passed.')
