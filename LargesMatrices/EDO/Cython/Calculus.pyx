import numpy as np
cimport numpy as np
from libc.math cimport sin, cos, pi, sqrt

cdef double beta = 0.1  # Set your value of beta here

cpdef double Sol(double x, double y):
    cdef double t2 = x ** 2
    cdef double t6 = cos(pi * y)
    return sin(pi * x + t6 * (t2 - x) * beta)

cpdef double[::1] champMagNormaliseRenverse(double t, double[::1] u):
    cdef double x = u[0]
    cdef double y = u[1]

    cdef double t1 = beta * pi
    cdef double t3 = -1 + x
    cdef double t4 = pi * y
    cdef double t5 = sin(t4)
    cdef double t7 = beta ** 2
    cdef double t8 = x ** 2
    cdef double t10 = pi * (t8 - x)
    cdef double t11 = 2 * x
    cdef double t15 = cos(t4)
    cdef double t16 = t15 ** 2
    cdef double t23 = pi ** 2
    cdef double t24 = t3 ** 2
    cdef double t30 = sqrt(-t16 * (t10 - t11 + 0.1) * (t10 + t11 - 0.1) * t7 + 0.4 * t15 * (x - 0.1 / 0.2) * t1 + (t7 * t24 * t8 + 0.1) * t23)
    cdef double MatBx = 0.1 / t30 * t5 * t3 * x * t1

    t1 = 2 * x
    t5 = cos(pi * y)
    t8 = beta ** 2
    t9 = x ** 2
    t11 = pi * (t9 - x)
    t15 = t5 ** 2
    t23 = pi ** 2
    t25 = (-1 + x) ** 2
    t31 = sqrt(-t15 * (t11 - t1 + 0.1) * (t11 + t1 - 0.1) * t8 + 0.4 * t5 * (x - 0.1 / 0.2) * pi * beta + (t8 * t25 * t9 + 0.1) * t23)
    cdef double MatBy = 0.1 / t31 * (t5 * (t1 - 1) * beta + pi)
    return np.array([MatBx, MatBy])
