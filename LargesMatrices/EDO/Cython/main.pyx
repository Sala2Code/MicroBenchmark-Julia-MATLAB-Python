import numpy as np
cimport numpy as np
import timeit
import cython

from globalVar import *
from Flots import *
from Calculus import *

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef main():
    cdef int NpointsX = 15
    cdef int NpointsY = 15

    cdef double[:] X = np.linspace(0, 1, NpointsX)
    cdef double[:] Y = np.linspace(0, 1, NpointsY)

    cdef double[:, :] solTheorique = np.zeros((NpointsY, NpointsX))
    cdef double[:, :] solPropagee = np.zeros((NpointsY, NpointsX))
    cdef double[:, :] solError = np.zeros((NpointsY, NpointsX))

    cdef Py_ssize_t i, j
    cdef double xIci, yIci
    cdef double[:] foot

    for i in range(NpointsX-2, 0, -1):
        for j in range(NpointsY-2, 0, -1):
            xIci = X[i]
            yIci = Y[j]
            # % changement i,j devient j,i
            solTheorique[j, i] = Sol(xIci, yIci)

            _, foot = propagation(np.array([xIci, yIci]))
            solPropagee[j, i] = Sol(foot[0], foot[1])

            solError[j, i] = abs(solPropagee[j, i] - solTheorique[j, i])


def timeMain():
    print(timeit.timeit(main, number=10)/10)