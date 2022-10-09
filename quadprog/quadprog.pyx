import numpy as np
from libc.stdio cimport printf

cdef extern nogil:
    int qpgen2_(double *G, double *av, int n,
                double *xv, double *lagr, double *obj,
                double *C, double *bv, int q, int meq,
                int* iact, int* nact, int* iter,
                double* work, int factorized)


def solve_qp(double[:, :] G, double[:] a, double[:, :] C=None, double[:] b=None, int meq=0, factorized=False):
    """Solve a strictly convex quadratic program

    Minimize     1/2 x^T G x - a^T x
    Subject to   C.T x >= b

    This routine uses the the Goldfarb/Idnani dual algorithm [1].

    References
    ---------
    ... [1] D. Goldfarb and A. Idnani (1983). A numerically stable dual
        method for solving strictly convex quadratic programs.
        Mathematical Programming, 27, 1-33.

    Parameters
    ----------
    G : array, shape=(n, n)
        matrix appearing in the quadratic function to be minimized
    a : array, shape=(n,)
        vector appearing in the quadratic function to be minimized
    C : array, shape=(n, m)
        matrix defining the constraints under which we want to minimize the
        quadratic function
    b : array, shape=(m), default=None
        vector defining the constraints
    meq : int, default=0
        the first meq constraints are treated as equality constraints,
        all further as inequality constraints (defaults to 0).
    factorized : bool, default=False
        If True, then we are passing :math:`R^{−1}` instead of the matrix G
        in the argument G, where :math:`G = R^T R` and R is upper triangular.

    Returns
    -------
    x : array, shape=(n,)
        vector containing the solution of the quadratic programming problem.
    f : float
        the value of the quadratic function at the solution.
    xu : array, shape=(n,)
        vector containing the unconstrained minimizer of the quadratic function
    iterations : tuple
        2-tuple. the first component contains the number of iterations the
        algorithm needed, the second indicates how often constraints became
        inactive after becoming active first.
    lagrangian : array, shape=(m,)
        vector with the Lagragian at the solution.
    iact : array
        vector with the indices of the active constraints at the solution.
    """

    cdef int n1, n2, m1
    n1, n2 = G.shape[0], G.shape[1]

    if C is None and b is None:
        C = np.zeros((n1, 1))
        b = -np.ones(1)
        meq = 0

    n3, m1 = C.shape[0], C.shape[1]
    if n1 != n2:
        raise ValueError('G must be a square matrix. Receive shape=(%d,%d)' % (n1, n2))
    if len(a) != n1:
        raise ValueError('G and a must have the same dimension. Received G as (%d,%d) and a as (%d,)' % (n1, n2, len(a)))
    if n1 != n3:
        raise ValueError('G and C must have the same first dimension. Received G as (%d,%d) and C as (%d, %d)' % (n1, n2, n3, m1))
    if len(b) != m1:
        raise ValueError('The number of columns of C must match the length of b. Received C as (%d, %d) and b as (%d,)' % (n2, m1, len(b)))

    cdef double[::1, :] G_ = np.array(G, copy=True, order='F')
    cdef double[::1, :] C_ = np.array(C, copy=True, order='F')
    cdef double[::1] a_ = np.array(a, copy=True, order='F')
    cdef double[::1] b_ = np.array(b, copy=True, order='F')
    cdef double[::1] sol = np.zeros(n1)
    cdef double[::1] lagr = np.zeros(m1)
    cdef double crval = 0
    cdef int meq_ = meq
    cdef int[::1] iact = np.zeros(m1, dtype=np.int32)
    cdef int nact = 0
    cdef int[::1] iters = np.zeros(2, dtype=np.int32)
    cdef double[::1] work = np.zeros(2*n1 + 2*m1 + min(n1, m1)*(min(n1, m1)+5)//2)
    cdef int factorized_ = 1 if factorized else 0

    cdef int result
    with nogil:
        result = qpgen2_(&G_[0, 0], &a_[0], n1,
                         &sol[0], &lagr[0], &crval,
                         &C_[0, 0], &b_[0], m1, meq_,
                         &iact[0], &nact, &iters[0],
                         &work[0], factorized_)

    if result == 1:
        raise ValueError('constraints are inconsistent, no solution')
    if result == 2:
        raise ValueError('matrix G is not positive definite')

    return np.asarray(sol), crval, np.asarray(a_), np.asarray(iters), np.asarray(lagr), np.asarray(iact)[:nact]
