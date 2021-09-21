import numpy as np
import scipy.optimize
import scipy.stats
from quadprog import solve_qp


def solve_qp_scipy(G, a, C, b, meq=0):
    # Minimize     1/2 x^T G x - a^T x
    # Subject to   C.T x >= b
    def f(x):
        return 0.5 * np.dot(x, G).dot(x) - np.dot(a, x)

    if C is not None and b is not None:
        constraints = [{
            'type': 'ineq',
            'fun': lambda x, C=C, b=b, i=i: (np.dot(C.T, x) - b)[i]
        } for i in range(C.shape[1])]
    else:
        constraints = []

    result = scipy.optimize.minimize(
        f, x0=np.zeros(len(G)), method='COBYLA', constraints=constraints,
        tol=1e-10, options={'maxiter': 2000})
    return result


def verify(G, a, C=None, b=None):
    xf, f, xu, iters, lagr, iact = solve_qp(G, a, C, b)
    result = solve_qp_scipy(G, a, C, b)
    np.testing.assert_array_almost_equal(result.x, xf)
    np.testing.assert_array_almost_equal(result.fun, f)


def test_1():
    G = np.eye(3, 3)
    a = np.array([0, 5, 0], dtype=np.double)
    C = np.array([[-4, 2, 0], [-3, 1, -2], [0, 0, 1]], dtype=np.double)
    b = np.array([-8, 2, 0], dtype=np.double)
    xf, f, xu, iters, lagr, iact = solve_qp(G, a, C, b)
    np.testing.assert_array_almost_equal(xf, [0.4761905, 1.0476190, 2.0952381])
    np.testing.assert_almost_equal(f, -2.380952380952381)
    np.testing.assert_almost_equal(xu, [0, 5, 0])
    np.testing.assert_array_equal(iters, [3, 0])
    np.testing.assert_array_almost_equal(lagr, [0.0000000, 0.2380952, 2.0952381])

    verify(G, a, C, b)


def test_2():
    G = np.eye(3, 3)
    a = np.array([0, 0, 0], dtype=np.double)
    C = np.ones((3, 1))
    b = -1000 * np.ones(1)
    verify(G, a, C, b)
    verify(G, a)


def test_3():
    random = np.random.RandomState(0)
    G = scipy.stats.wishart(scale=np.eye(3, 3), seed=random).rvs()
    a = random.randn(3)
    C = random.randn(3, 2)
    b = random.randn(2)
    verify(G, a, C, b)
    verify(G, a)


def test_4():
    n = 66
    X = np.full((n, n), 1e-20)
    X[np.diag_indices_from(X)] = 1.0
    y = np.arange(n, dtype=float)

    random = np.random.RandomState(1)
    G = np.dot(X.T, X)
    a = np.dot(X, y)
    C = np.identity(n)
    b = y + random.rand(n)
    
    verify(G, a, C, b)

