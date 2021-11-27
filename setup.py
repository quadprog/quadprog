from setuptools import setup, Extension

long_description = """Minimize     1/2 x^T G x - a^T x

Subject to   C.T x >= b

This routine uses the the Goldfarb/Idnani dual algorithm [1].

References
----------
1) D. Goldfarb and A. Idnani (1983). A numerically stable dual
   method for solving strictly convex quadratic programs.
   Mathematical Programming, 27, 1-33.
"""

classifiers = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Education",
    "Intended Audience :: Financial and Insurance Industry",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Mathematics"
]

extensions = [Extension('quadprog', [
    'quadprog/linear-algebra.c',
    'quadprog/qr-update.c',
    'quadprog/quadprog.c',
    'quadprog/solve.QP.c',
])]

setup(
    name="quadprog",
    version="0.1.11",
    description="Quadratic Programming Solver",
    long_description=long_description,
    url="https://github.com/quadprog/quadprog",
    author="Robert T. McGibbon",
    author_email='rmcgibbo@gmail.com',
    license='GPLv2+',
    classifiers=classifiers,
    ext_modules=extensions,
    install_requires=["numpy"]
)
