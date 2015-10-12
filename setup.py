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

from setuptools import setup, Extension
from Cython.Build import cythonize

##########################
VERSION = "0.1.6"
__version__ = VERSION
##########################


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
                           'quadprog/aind.c', 'quadprog/solve.QP.c',
                           'quadprog/util.c', 'quadprog/dpofa.c',
                           'quadprog/daxpy.c', 'quadprog/ddot.c',
                           'quadprog/dscal.c', 'quadprog/f2c_lite.c'],
             include_dirs=['quadprog'], language='c++')
]

setup(
    name='quadprog',
    author="Robert T. McGibbon",
    author_email='rmcgibbo@gmail.com',
    url="https://github.com/rmcgibbo/quadprog",
    description=DOCLINES[0],
    version=__version__,
    long_description="\n".join(DOCLINES[2:]),
    license='GPLv2+',
    zip_safe=False,
    ext_modules=cythonize(extensions),
)
