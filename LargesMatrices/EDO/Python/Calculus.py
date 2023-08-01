import numpy as np
from globalVar import *

def Sol(x,y):
    t2 = x ** 2
    t6 = np.cos(np.pi * y)
    return np.sin(np.pi * x + t6 * (t2 - x) * beta)


def champMagNormaliseRenverse(t, u):
    x = u[0]
    y = u[1]
    
    t1 = beta * np.pi
    t3 = -1 + x
    t4 = np.pi * y
    t5 = np.sin(t4)
    t7 = beta ** 2
    t8 = x ** 2
    t10 = np.pi * (t8 - x)
    t11 = 2 * x
    t15 = np.cos(t4)
    t16 = t15 ** 2
    t23 = np.pi ** 2
    t24 = t3 ** 2
    t30 = np.sqrt(-t16 * (t10 - t11 + 0.1e1) * (t10 + t11 - 0.1e1) * t7 + 0.4e1 * t15 * (x - 0.1e1 / 0.2e1) * t1 + (t7 * t24 * t8 + 0.1e1) * t23)
    MatBx = 0.1e1 / t30 * t5 * t3 * x * t1
    
    t1 = 2 * x
    t5 = np.cos(np.pi * y)
    t8 = beta ** 2
    t9 = x ** 2
    t11 = np.pi * (t9 - x)
    t15 = t5 ** 2
    t23 = np.pi ** 2
    t25 = (-1 + x) ** 2
    t31 = np.sqrt(-t15 * (t11 - t1 + 0.1e1) * (t11 + t1 - 0.1e1) * t8 + 0.4e1 * t5 * (x - 0.1e1 / 0.2e1) * np.pi * beta + (t8 * t25 * t9 + 0.1e1) * t23)
    MatBy = 0.1e1 / t31 * (t5 * (t1 - 1) * beta + np.pi)
    
    # b = [aBx(x,y) ; aBy(x,y)];
    # b = [MatBx ; MatBy];
    
    return np.reshape(np.array([[-MatBx], [-MatBy]]) , (2,))
