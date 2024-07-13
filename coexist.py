import numpy as np
import scipy.optimize, scipy.spatial, scipy.special
import math
from free_energy import fe

## Calculate the m coexisting phases for a N-component system
## m(<=N) can be greater than the number of coexistence phases. If that is the case, set conver=False for the first round of optimization. 
## Then merge the overlapping phase points and do a second round with conver=True.

def coex_exact(phiguess, muguess, cguess, target, verbose=True, conver=False):
    tol = 10**(-12)
    ## phiguess: initial guess for concentration of the m coexisting phases (named 1 to m), array size of [m,N]
    ## muguess: initial guess for the chemical potential at coexistence, array size of N
    ## cguess: initial guess for the volume fractions of phases (2 to m), array size of m-1. The volume fraction of phase 1 is 1-np.sum(cguess)
    ## target: target concentration for phase-coexistence calculation, array size of N
    ## tot: precision of calculation

    def diff(muc,target,calcphi=False):
        ## This is the multivariable function that the solver tries to minimize.
        ## muc: combined array of chemical potential (first N elements) and volume fraction, array size of N+m-1
        mu = muc[:len(muguess)]
        volf = muc[len(muguess):]

        ## Grand potential as a function of phi vector (size N)
        def Omega(phi):
            if np.any(phi < 0):
                return np.inf
            return fe(phi)-np.dot(mu,phi)

        res0 = scipy.optimize.minimize(Omega,
                                       phiguess[0], method='Nelder-Mead',
                                       options={'xatol': tol, 'fatol': tol})
        ## Find the minimum of grand potential near the first guessed phi vector, and the corresponding grand potential value
        phi0 = res0.x
        Omega0 = res0.fun
        phimin, Omegamin = [], []
        ## Similarly, find the minima near the other initial phi vectors.
        for phi in phiguess[1:]:
            res = scipy.optimize.minimize(Omega,
                                          phi, method='Nelder-Mead',
                                          options={'xatol': tol, 'fatol': tol})
            phimin.append(res.x)
            Omegamin.append(res.fun)
        ## The returned vector d also has the size of m+N-1
        ## The first m-1 elements record the grand potential differences referenced by the first phase.
        d = np.zeros(len(muguess)+len(cguess))
        d[:len(phimin)] = np.array(Omegamin) - Omega0
        d[len(phimin):] = phiconstraints(volf,phi0,phimin,target)
        ## The last N elements are the residue errors in concentration.

        ## Set calcphi=True to record coexisting phi after optimization.
        if calcphi:
            print("  --> |diff| = %10g, diff =" % np.linalg.norm(d), d)
            return [phi0] + phimin
        return d

    ## Combine arrays to build the initial guess vector 
    mucguess = np.concatenate((muguess,cguess))

    ## solve for phase coexistence using the modified Powell method
    res = scipy.optimize.root(diff, mucguess, args=(target,),method='hybr')
    if verbose: print("COEX STATUS:", res.success, res.message)

    ## If the results are unconverged but we want to get something, try Levenberg–Marquardt method.
    if conver and res.success==False:
        print('using LM')
        res = scipy.optimize.root(diff, mucguess, args=(target,),method='lm')
        if verbose: print("COEX STATUS:", res.success, res.message)
    mu = res.x[:len(muguess)]
    volf = res.x[len(muguess):]
    muc = np.concatenate((mu,volf))
    phi = diff(muc, target, calcphi=True)
    ## phi: coexisting phases, array size of [m,N]
    ## mu: chemical potential, array size of N
    ## volf: volume fractions of phases (NO.2 to NO.m), array size of N-1
    return phi, mu, volf

## Calculate the residue errors in concentration
def phiconstraints(volf,phi0,phimin,target):
    resi = phi0*(1-np.sum(volf)) + np.dot(volf,np.array(phimin)) - target
    return resi

