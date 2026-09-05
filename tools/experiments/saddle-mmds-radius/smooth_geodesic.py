"""Smooth saddle geodesics by adaptive shooting, with analytic sensitivities.

Coordinates are in the unit base disk; the surface is (r*x,r*y,r²*(x²-y²)).
Negative Gaussian curvature, completeness and simple connectivity make the
converged geodesic unique and globally minimizing. See README.md.
"""
import numpy as np
from numba import njit, prange


@njit(cache=True)
def rhs(y,k):
    x,z,v,w=y[:4]
    den=1+k*(x*x+z*z)
    a=k*(v*v-w*w)/den
    ax=-2*k*x*a/den;az=-2*k*z*a/den
    av=2*k*v/den;aw=-2*k*w/den
    out=np.empty(12)
    out[0]=v;out[1]=w;out[2]=-a*x;out[3]=a*z
    for j in range(2):
        b=4+4*j
        dx,dz,dv,dw=y[b:b+4]
        da=ax*dx+az*dz+av*dv+aw*dw
        out[b]=dv;out[b+1]=dw
        out[b+2]=-a*dx-x*da;out[b+3]=a*dz+z*da
    return out


@njit(cache=True)
def integrate(start,velocity,r,rtol=2e-10,atol=2e-12):
    """Dormand-Prince 5(4), integrating position/velocity plus 2 sensitivities."""
    y=np.zeros(12);y[:2]=start;y[2:4]=velocity;y[6]=1;y[11]=1
    t=0.;step=.02;k=4*r*r;steps=0
    while t<1 and steps<20000:
        step=min(step,1-t)
        k1=rhs(y,k)
        k2=rhs(y+step*k1/5,k)
        k3=rhs(y+step*(3*k1/40+9*k2/40),k)
        k4=rhs(y+step*(44*k1/45-56*k2/15+32*k3/9),k)
        k5=rhs(y+step*(19372*k1/6561-25360*k2/2187+64448*k3/6561-212*k4/729),k)
        k6=rhs(y+step*(9017*k1/3168-355*k2/33+46732*k3/5247+49*k4/176-5103*k5/18656),k)
        trial=y+step*(35*k1/384+500*k3/1113+125*k4/192-2187*k5/6784+11*k6/84)
        k7=rhs(trial,k)
        err=step*(71*k1/57600-71*k3/16695+71*k4/1920-17253*k5/339200+22*k6/525-k7/40)
        error=0.
        for j in range(12):
            scale=atol+rtol*max(abs(y[j]),abs(trial[j]))
            error=max(error,abs(err[j])/scale)
        if not np.isfinite(error):return y,False,steps
        if error<=1:
            y=trial;t+=step
        factor=5. if error==0 else min(5.,max(.1,.9*error**(-.2)))
        step*=factor;steps+=1
        if step<1e-14:return y,False,steps
    return y,t>=1,steps


@njit(cache=True)
def shoot(a,b,r,initial,tolerance=2e-9,rtol=2e-10):
    velocity=initial.copy();total_steps=0
    for iteration in range(35):
        y,ok,steps=integrate(a,velocity,r,rtol,rtol/100);total_steps+=steps
        if not ok:return velocity,np.inf,False,iteration,total_steps
        residual=y[:2]-b;error=np.sqrt(np.sum(residual*residual))
        if error<tolerance:
            length=r*np.sqrt(np.sum(velocity**2)+4*r*r*(a[0]*velocity[0]-a[1]*velocity[1])**2)
            return velocity,length,True,iteration,total_steps
        j00=y[4];j10=y[5];j01=y[8];j11=y[9]
        det=j00*j11-j01*j10
        if abs(det)<1e-16:return velocity,np.inf,False,iteration,total_steps
        delta=np.array([(j11*residual[0]-j01*residual[1])/det,(-j10*residual[0]+j00*residual[1])/det])
        accepted=False;scale=1.
        for line in range(16):
            trial=velocity-scale*delta
            check,ok,steps=integrate(a,trial,r,rtol,rtol/100);total_steps+=steps
            if ok and np.linalg.norm(check[:2]-b)<error:
                velocity=trial;accepted=True;break
            scale*=.5
        if not accepted:return velocity,np.inf,False,iteration,total_steps
    return velocity,np.inf,False,35,total_steps


@njit(cache=True,parallel=True)
def shoot_batch(a,b,r,initial,rtol=2e-10):
    n=len(a);velocities=np.empty_like(initial);lengths=np.empty(n)
    success=np.zeros(n,np.bool_);iterations=np.zeros(n,np.int64);steps=np.zeros(n,np.int64)
    for j in prange(n):
        v,d,ok,it,st=shoot(a[j],b[j],r,initial[j],2e-9,rtol)
        velocities[j]=v;lengths[j]=d;success[j]=ok;iterations[j]=it;steps[j]=st
    return velocities,lengths,success,iterations,steps


def pair_coordinates(u):
    i,j=np.triu_indices(len(u),1)
    return u[i],u[j]


def solve_pairs(a,b,r,initial,previous=None,rtol=2e-10):
    """Retry failed Newton solves by finer continuation from the previous case."""
    v,d,ok,it,steps=shoot_batch(a,b,r,initial,rtol)
    retries=0
    if not np.all(ok):
        if previous is None:
            old_a,old_b,old_r=a,b,.01
            old_v,_,old_ok,_,_=shoot_batch(a,b,.01,b-a,rtol)
            assert np.all(old_ok)
        else:old_a,old_b,old_r,old_v=previous
        for divisions in [2,4,8,16,32,64]:
            indices=np.where(~ok)[0]
            if len(indices)==0:break
            warm=old_v[indices].copy();working=np.ones(len(indices),bool)
            for j in range(1,divisions+1):
                fraction=j/divisions
                aa=(1-fraction)*old_a[indices]+fraction*a[indices]
                bb=(1-fraction)*old_b[indices]+fraction*b[indices]
                rr=np.exp((1-fraction)*np.log(old_r)+fraction*np.log(r))
                warm,dd,good,ii,ss=shoot_batch(aa,bb,rr,warm,rtol)
                working &= good;steps[indices]+=ss;retries+=len(indices)
                if not np.any(working):break
            fixed=indices[working]
            v[fixed]=warm[working];d[fixed]=dd[working];ok[fixed]=True
    if not np.all(ok):raise RuntimeError(f'Unresolved smooth geodesics: {np.where(~ok)[0].tolist()} at r={r}')
    return v,d,dict(retry_solves=retries,max_newton_iterations=int(it.max()),rk_steps=int(steps.sum()))


def validate():
    """Validate against independent SciPy DOP853 integration and collocation."""
    from scipy.integrate import solve_ivp, solve_bvp, quad
    rng=np.random.default_rng(7201)
    angles=rng.uniform(0,2*np.pi,24);rad=np.sqrt(rng.uniform(.02,.95,24))
    points=rad[:,None]*np.c_[np.cos(angles),np.sin(angles)]
    a=points[:12];b=points[12:];warm=b-a;previous=None;rows=[]
    for r in [.1,.25,.5,1,2,4,8,16,32,64]:
        velocity,length,_=solve_pairs(a,b,r,warm,previous)
        if r in [.1,1,4,16,64]:
            for j in range(len(a)):
                k=4*r*r
                def ode(t,state):
                    x,y,v,w=state
                    accel=k*(v*v-w*w)/(1+k*(x*x+y*y))
                    return np.array([v,w,-accel*x,accel*y])
                ivp=solve_ivp(ode,[0,1],np.r_[a[j],velocity[j]],method='DOP853',rtol=2e-12,atol=2e-13,dense_output=True)
                endpoint=np.linalg.norm(ivp.y[:2,-1]-b[j])
                mesh=np.unique(np.r_[ivp.t,np.linspace(0,1,60)])
                bvp=solve_bvp(ode,lambda ya,yb:np.r_[ya[:2]-a[j],yb[:2]-b[j]],mesh,ivp.sol(mesh),tol=1e-6,max_nodes=30000)
                assert bvp.success, (r,j,bvp.message,len(bvp.x))
                v0=bvp.sol(0)[2:]
                exact=r*np.sqrt(np.dot(v0,v0)+k*(a[j,0]*v0[0]-a[j,1]*v0[1])**2)
                state=ivp.sol(mesh);energy=np.sum(state[2:]**2,axis=0)+k*(state[0]*state[2]-state[1]*state[3])**2
                rows.append(dict(radius=r,pair=j,endpoint_error_over_radius=endpoint,
                                 relative_distance_error_bvp=abs(exact-length[j])/length[j],
                                 relative_energy_drift=float(np.ptp(energy)/energy[0])))
        previous=(a,b,r,velocity);warm=velocity
    for r in [.1,1,8,64]:
        # A ruling y=x is a straight line in 3D and therefore exactly geodesic.
        a0=np.array([.6,.6]);b0=np.array([-.3,-.3])
        _,d,ok,_,_=shoot(a0,b0,r,b0-a0)
        assert ok and abs(d-r*np.linalg.norm(b0-a0))<1e-8
        # An axis meridian is fixed by a reflection; check its exact arc length.
        a0=np.array([.7,0]);b0=np.array([-.4,0])
        v,d,_=solve_pairs(a0[None,:],b0[None,:],r,(b0-a0)[None,:])
        exact=quad(lambda x:np.sqrt(1+4*x*x),-.4*r,.7*r,epsabs=1e-10)[0]
        assert abs(d[0]-exact)/exact<1e-8
    assert max(a['endpoint_error_over_radius'] for a in rows)<2e-5
    assert max(a['relative_distance_error_bvp'] for a in rows)<1e-7
    return rows
