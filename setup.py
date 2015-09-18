"""Solve a Quadratic Programming Problem

Authors
-------
-  Robert T. McGibbon rmcgibbo@gmail.com
"""
import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import numpy as np
from Cython.Build import cythonize
from numpy.distutils.command import build_src
class build_src(build_src.build_src):
    def f2py_sources(self, sources, extension):
        return sources


from numpy.distutils.system_info import get_info
blas_opt = get_info('blas_opt', notfound_action=2)

# libs = [blas_opt['libraries'][0], lapack_opt['libraries'][0]]
# lib_dirs = [blas_opt['library_dirs'][0], lapack_opt['library_dirs'][0]]
#
# print(libs, lib_dirs)

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
    Extension('quadprog', ['quadprog/quadprog.pyx',
                           'quadprog/aind.f', 'quadprog/solve.QP.f',
                           'quadprog/util.f', 'quadprog/dpofa.f',
                           'quadprog/fortranwrapper.c'],
              include_dirs=[np.get_include(), 'quadprog'],
              define_macros=blas_opt['define_macros'],
              extra_compile_args=blas_opt['extra_compile_args'],
              extra_link_args=blas_opt['extra_link_args'],
              language='c++'),
]

setup(
    name='quadprog',
    author="Robert T. McGibbon",
    author_email='rmcgibbo@gmail.com',
    cmdclass={'build_src': build_src},
    url="https://github.com/rmcgibbo/quadprog",
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    license='GPL',
    install_requires=['numpy'],
    zip_safe=False,
    ext_modules=cythonize(extensions),
)
