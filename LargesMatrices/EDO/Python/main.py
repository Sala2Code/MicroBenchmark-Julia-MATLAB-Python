import numpy as np
import timeit

from globalVar import *
from Flots import *
from Calculus import *

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

            _, foot = propagation([xIci, yIci])
            solPropagee[j, i] = Sol(foot[0], foot[1])

            solError[j,i] = abs(solPropagee[j,i] - solTheorique[j,i])


print(timeit.timeit(main, number=10)/10)

# main()