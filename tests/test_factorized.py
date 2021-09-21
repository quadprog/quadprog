import numpy as np
import scipy.stats
import scipy.linalg
from quadprog import solve_qp
random = np.random.RandomState(0)


def verify(G, a, C=None, b=None):
    xf0, f0 = solve_qp(G, a, C, b)[0:2]
    xf1, f1 = solve_qp(scipy.linalg.inv(scipy.linalg.cholesky(G)), a, C, b, factorized=True)[0:2]

    np.testing.assert_array_almost_equal(xf0, xf1)
    np.testing.assert_almost_equal(f0, f1)


def test_1():
    G = scipy.stats.wishart(scale=np.eye(3,3), seed=random).rvs()
    a = random.randn(3)
    verify(G, a)


def test_2():
    G = scipy.stats.wishart(scale=np.eye(3,3), seed=random).rvs()
    a = random.randn(3)
    C = random.randn(3, 2)
    b = random.randn(2)
    verify(G, a, C, b)
