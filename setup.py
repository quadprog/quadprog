"""Solve a Quadratic Programming Problem

Authors
-------
-  Robert T. McGibbon rmcgibbo@gmail.com
"""
from setuptools import find_packages, setup, Extension
import numpy as np
from Cython.Distutils import build_ext

DOCLINES = __doc__.split("\n")
CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GPL
Programming Language :: Python
Operating System :: OS Independent
"""

extensions = [
    Extension('quadprog', ['quadprog/quadprog.pyx'],
              include_dirs=[np.get_include()],
              language='c++'),
]

setup(
    name='quadprog',
    author="Robert T. McGibbon",
    author_email='rmcgibbo@gmail.com',
    cmdclass={'build_ext': build_ext},
    url="https://github.com/rmcgibbo/quadprog",
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    license='GPL',
    install_requires=['numpy'],
    zip_safe=False,
    ext_modules=extensions,
)
