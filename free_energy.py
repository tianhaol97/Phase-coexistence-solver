import numpy as np
import math


### This file provides the setup of an example of phase-coexistence calculation
### Number of explicit component is 3
nc = 3

kai= np.array([[0.,9.,3.,4.5],
               [9.,0.,3.,4.5],
               [3.,3.,0.,1.5],
               [4.5,4.5,1.5,0.]])


def fe(rho):
    rhoe = np.zeros(nc+1)
    rhoe[:nc] = np.copy(rho)
    rhoe[-1] = 1-np.sum(rho)
    if np.all(rhoe>=-10**(-9)):
        f = np.dot(np.matmul(rhoe,kai),rhoe)/2
        for i in range(nc+1):
            if rhoe[i]>0:
                f = f+rhoe[i]*math.log(rhoe[i])
    else:
        f = np.inf
    return f


