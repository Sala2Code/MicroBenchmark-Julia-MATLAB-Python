from scipy.integrate import solve_ivp
import numpy as np
cimport numpy as np
import cython

from Calculus import *
from globalVar import *

ctypedef np.float64_t[:] arr_t # définir hors de la fonction 

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double event(double t, arr_t X):
    cdef double event_value = X[1] * (X[1] - 1)
    event_value = 0 if event_value < 0 else event_value
    return event_value

def propagation(Ydepart): # possible optimisation mais je n'ai pas réussi, assez difficile à faire.
    tmax = 15.0
    intervalleT = (0.0, tmax)
    solution = solve_ivp(champMagNormaliseRenverse, intervalleT, Ydepart, method='RK45', events=event, rtol=rTol, atol=aTol)
    return solution.t_events[-1][0], solution.y_events[0][-1].T
