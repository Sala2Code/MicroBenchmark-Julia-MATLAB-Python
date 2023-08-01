from scipy.integrate import solve_ivp
import numpy as np

from Calculus import *
from globalVar import *

def event(t, X):
    event_value = X[1]*(X[1]-1)
    event_value = 0 if event_value < 0 else event_value
    return event_value

def propagation(Ydepart):
    tmax = 15.0
    intervalleT = (0.0, tmax)
    solution = solve_ivp(champMagNormaliseRenverse, intervalleT, Ydepart, method='RK45', events=event, rtol=rTol, atol=aTol)
    return solution.t_events[-1][0], solution.y_events[0][-1].T

