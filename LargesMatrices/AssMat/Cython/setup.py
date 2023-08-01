from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np

sourcefiles = ['ParamBeta.pyx', "Calculs.pyx","Lap.pyx", "lexico.pyx"]

setup(
    name='PlasmaCython',
    ext_modules=cythonize(sourcefiles),
    zip_safe=False,
    include_dirs=[np.get_include()]
)


# * To run
# python setup.py build_ext --inplace
# python main.py