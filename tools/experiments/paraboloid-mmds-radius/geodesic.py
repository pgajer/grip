"""Smooth geodesic distances on z=x²+y², without a sampled graph.

The metric is (1+4 rho²) d rho² + rho² d theta². Clairaut's
constant c is the minimum radius of a turning geodesic. Closed-form
angle/length primitives reduce the inverse problem to scalar bisection.
See GEODESIC.md for branch selection and the shortest-path argument.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.spatial.distance import squareform


def primitives(c, t):
    """Angle and length from the turn to rho=sqrt(c²+t²)."""
    a = 1 + 4*c*c
    h = np.sqrt(a+4*t*t)
    v = np.arcsinh(2*t/np.sqrt(a))
    angle = 2*c*v + np.arctan2(t, c*h)
    length = .5*t*h + .25*a*v
    return angle, length


def pair_distances(a, b, angle, iterations=62, details=False):
    """Shortest distances for broadcastable radii and angular gaps in [0,pi]."""
    a,b,angle = np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float),np.asarray(angle,float))
    a,b = np.minimum(a,b),np.maximum(a,b)
    gap = np.remainder(np.abs(angle),2*np.pi)
    gap = np.minimum(gap,2*np.pi-gap)
    # With c=a sin(phi), t_a=a cos(phi), avoid subtracting nearly equal
    # squares for short turning paths near a parallel.
    tb0 = np.sqrt(np.maximum((b-a)*(b+a),0))
    threshold = primitives(a,tb0)[0]
    turning = gap > threshold
    lo = np.zeros(a.shape); hi = np.full(a.shape,np.pi/2)
    for _ in range(iterations):
        phi = .5*(lo+hi)
        c=a*np.sin(phi); ta=a*np.cos(phi)
        tb=np.sqrt(np.maximum((b-a)*(b+a)+ta*ta,0))
        aa, _=primitives(c,ta); ab, _=primitives(c,tb)
        predicted=np.where(turning,ab+aa,ab-aa)
        move_lo=np.where(turning,predicted>gap,predicted<gap)
        lo=np.where(move_lo,phi,lo); hi=np.where(move_lo,hi,phi)
    phi=.5*(lo+hi); c=a*np.sin(phi); ta=a*np.cos(phi)
    tb=np.sqrt(np.maximum((b-a)*(b+a)+ta*ta,0))
    aa,la=primitives(c,ta); ab,lb=primitives(c,tb)
    length=np.where(turning,lb+la,lb-la)
    predicted=np.where(turning,ab+aa,ab-aa)
    # Meridian and pole distances are known exactly.
    radial = primitives(np.zeros_like(a),b)[1]-primitives(np.zeros_like(a),a)[1]
    length=np.where((gap==0)|(a==0),radial,length)
    length=np.maximum(length,0)
    if details:
        return length,dict(c=c,ta=ta,tb=tb,turning=turning,predicted_angle=predicted,
                           target_angle=gap,a=a,b=b)
    return length


def distance_matrix(rho, theta):
    i,j=np.triu_indices(len(rho),1)
    return squareform(pair_distances(rho[i],rho[j],theta[i]-theta[j]))


def validate():
    """Independent quadrature and Cartesian geodesic-ODE endpoint checks."""
    rng=np.random.default_rng(4701)
    records=[]
    for r in [.1,1.,4.,16.,64.]:
        for j in range(10):
            a,b=np.sort(r*rng.uniform(.08,1.,2))
            angle=rng.uniform(.02,np.pi-.02)
            length,info=pair_distances(a,b,angle,details=True)
            c=float(info['c']); ta=float(info['ta']); tb=float(info['tb'])
            turn=bool(info['turning'])
            intervals=[(0,ta),(0,tb)] if turn else [(ta,tb)]
            numeric_length=sum(quad(lambda t:np.sqrt(1+4*c*c+4*t*t),l,h,epsabs=1e-10,epsrel=1e-11)[0] for l,h in intervals)
            # An independent angular integral in phi=atan(t/c).
            numeric_angle=sum(quad(lambda p:np.sqrt(1+4*c*c/(np.cos(p)**2)),np.arctan2(l,c),np.arctan2(h,c),epsabs=1e-11,epsrel=1e-11)[0] for l,h in intervals)
            # Cartesian equation for the graph metric: x''=-4 x |x'|²/(1+4|x|²).
            sign=-1 if turn else 1
            initial=[a,0,sign*(ta/a)/np.sqrt(1+4*a*a),c/a]
            def ode(_,state):
                x=state[:2]; v=state[2:]
                return np.r_[v,-4*x*np.dot(v,v)/(1+4*np.dot(x,x))]
            sol=solve_ivp(ode,[0,float(length)],initial,rtol=3e-11,atol=1e-12)
            expected=np.array([b*np.cos(angle),b*np.sin(angle)])
            endpoint_error=np.linalg.norm(sol.y[:2,-1]-expected)/max(r,1e-12)
            records.append(dict(radius=r,case=j,turning=turn,length=float(length),
                                length_quad_relative_error=abs(numeric_length-length)/max(float(length),1e-12),
                                angle_quad_error=abs(numeric_angle-angle),
                                ode_endpoint_error_over_radius=float(endpoint_error)))
    # Equal endpoints, apex, same meridian, and opposite meridians.
    assert float(pair_distances(3,3,0)) == 0
    assert abs(float(pair_distances(0,3,1.1))-float(primitives(0,3)[1])) < 1e-12
    assert abs(float(pair_distances(.3,3,0))-float(primitives(0,3)[1]-primitives(0,.3)[1])) < 1e-12
    assert abs(float(pair_distances(.1,.1,np.pi))-2*float(primitives(0,.1)[1])) < 1e-10
    assert float(pair_distances(10,10,np.pi)) < 2*float(primitives(0,10)[1])
    assert max(a['length_quad_relative_error'] for a in records)<1e-10
    assert max(a['angle_quad_error'] for a in records)<1e-8
    assert max(a['ode_endpoint_error_over_radius'] for a in records)<1e-7
    return records
