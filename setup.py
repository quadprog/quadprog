"""Quadratic Programming Solver

Minimize     1/2 x^T G x - a^T x

Subject to   C.T x >= b

This routine uses the the Goldfarb/Idnani dual algorithm [1].

References
----------
1) D. Goldfarb and A. Idnani (1983). A numerically stable dual
   method for solving strictly convex quadratic programs.
   Mathematical Programming, 27, 1-33.
"""
import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import numpy as np
from Cython.Build import cythonize
from numpy.distutils.command import build_src

##########################
VERSION = "0.1.4"
__version__ = VERSION
##########################


class build_src(build_src.build_src):
    def f2py_sources(self, sources, extension):
        return sources


DOCLINES = __doc__.split("\n")
CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)
Programming Language :: Python
Operating System :: OS Independent
"""


extensions = [
    Extension('quadprog', ['quadprog/quadprog.pyx',
                           'quadprog/aind.f', 'quadprog/solve.QP.f',
                           'quadprog/util.f', 'quadprog/dpofa.f',
                           'quadprog/daxpy.f', 'quadprog/ddot.f',
                           'quadprog/dscal.f',
                           'quadprog/fortranwrapper.c'],
              language='c++'),
]

setup(
    name='quadprog',
    author="Robert T. McGibbon",
    author_email='rmcgibbo@gmail.com',
    cmdclass={'build_src': build_src},
    url="https://github.com/rmcgibbo/quadprog",
    description=DOCLINES[0],
    version=__version__,
    long_description="\n".join(DOCLINES[2:]),
    license='GPLv2+',
    zip_safe=False,
    ext_modules=cythonize(extensions),
)
