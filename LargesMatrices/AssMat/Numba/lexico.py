from numba import jit

import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
warnings.filterwarnings("ignore", message=".*The keyword argument 'nopython=False'.*")


@jit(nopython=True)
def lexico(i, j, Nxy):
    return i*Nxy[1]+j
