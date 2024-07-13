import numpy as np
import scipy.optimize, scipy.spatial, scipy.special
import math
import networkx as nx
from free_energy import fe
from coexist import coex_exact,phiconstraints

### The example is from regular solution model (Li and Jacobs, 2024)
### Binding energies (in unit of kT) are: N-N:-1.5, N-B:-0.5, B-B:-0.5

### parent concentration, a vector representing volume fractions of N1, N2, B successively
target = [0.25,0.25,0.05]

### initial guess of 4 phase vectors
guess_phi = np.array([[1e-8,0.96610169,0.01694915], 
                [0.01694915,0.01694915,0.08474576], 
                [0.01694915,0.01694915,0.06779661], 
                [0.96610169,1e-8,0.01694915]])

guess_mu = np.array([0.0700487,0.0700487,-1.28929441])
guess_vf = np.array([0.45,0.05,0.25])


print('First round of coexistence calculation starts:')
phi, mu, vol = coex_exact(guess_phi,guess_mu,guess_vf,target)
print('*************')
print('Intermediate results:')
print('phi:',phi)
print('mu:',mu)
print('vol:',vol)
print('*************')

### Merge vectors within tolerance. 
### Here we use a graph-based method. Merging can be done with other methods or manually.

def combinephi(phi):
    tol = 10**(-6)
    G = nx.Graph()
    phi = np.array(phi)
    for i in range(4):
        G.add_node(i)
        for j in range(i,4):
            if np.linalg.norm(phi[i]-phi[j])<tol:
                G.add_edge(i, j)
    return G

G = combinephi(phi)
phis = []
vol = np.insert(vol,0,1-np.sum(vol))
vf = []
for subg in nx.connected_components(G):
    ln = list(subg)
    phis.append(phi[ln[0]])
    vf.append(sum(vol[i] for i in ln))

phis = np.array(phis)
vf = np.array(vf[1:])

### 3 (out of 4) phase vectors remain:
### phis: [[1.55128666e-04 9.71529808e-01 1.59903308e-02]
### [1.56982780e-02 1.56982774e-02 8.37873373e-02]
### [9.71529805e-01 1.55128641e-04 1.59903327e-02]]
### vf: [0.50163852 0.24918076]



print('Merging the same phase vectors')
print('phi:',phis)
print('vol:',vf)
print('*************')
print('Second round of coexistence calculation starts:')
phif, muf, volf = coex_exact(phis,mu,vf,target,conver=True)
print('*************')
print('Final results:')
print('phi:',phif)
print('mu:',muf)
print('vol:',volf)