from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np

sourcefiles = ["main.pyx", "Calculus.pyx", "Flots.pyx", "globalVar.pyx"]

setup(
    name='EDO',
    ext_modules=cythonize(sourcefiles),
    zip_safe=False,
    include_dirs=[np.get_include()]
)

# * To run
# python setup.py build_ext --inplace
# python main.py