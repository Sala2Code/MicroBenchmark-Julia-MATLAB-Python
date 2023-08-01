import numpy as np
import time
from numba import jit, objmode, float64

from globalVar import *
from Flots import *
from Calculus import *

@jit(nopython=True)
def main():
    NpointsX = 15
    NpointsY = 15

    X = np.linspace(0, 1, NpointsX)
    Y = np.linspace(0, 1, NpointsY)

    solTheorique = np.zeros((NpointsY, NpointsX))
    solPropagee = np.zeros((NpointsY, NpointsX))
    solError = np.zeros((NpointsY, NpointsX))

    for i in range(NpointsX-2,0,-1):
        for j in range(NpointsX-2,0,-1):
            xIci = X[i]
            yIci = Y[j]
            # % changement i,j devient j,i
            solTheorique[j,i] = Sol(xIci, yIci)
            with objmode(): # error with propagation...
                _, foot = propagation(np.array([xIci, yIci]))
                solPropagee[j, i] = Sol(foot[0], foot[1])

            solError[j,i] = abs(solPropagee[j,i] - solTheorique[j,i])

with objmode():
    t = 0
    n=1
    for i in range(n):
        t1 = time.time()
        main()  
        t+=time.time()-t1
    print(t/n)