import ctypes
import sys

import pytest

import quadprog


@pytest.mark.skipif(sys.platform == "win32", reason="nothing is exported on Windows")
@pytest.mark.parametrize("name", ["dot", "axpy", "scal", "cholesky", "qr_insert", "qpgen2_"])
def test_internal_symbols_not_exported(name):
    lib = ctypes.CDLL(quadprog.__file__)
    with pytest.raises(AttributeError):
        getattr(lib, name)
