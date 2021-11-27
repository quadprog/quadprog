/*
 * Solve a strictly convex quadratic program:
 *
 *  minimize     1/2 x^T G x - a^T x
 *  subject to   C1^T x  = b1
 *               C2^T x >= b2

 *  This routine uses the the Goldfarb/Idnani dual algorithm [1].

 *  References
 *  ---------
 *  ... [1] D. Goldfarb and A. Idnani (1983). A numerically stable dual
 *      method for solving strictly convex quadratic programs.
 *      Mathematical Programming, 27, 1-33.
 */

double fabs(double);
double sqrt(double);

void axpy(int n, double a, double x[], double y[]);
double dot(int n, double x[], double y[]);
void triangular_multiply(int n, double a[], double b[]);
void triangular_multiply_transpose(int n, double a[], double b[]);
void triangular_solve(int n, double a[], double b[]);
void triangular_solve_transpose(int n, double a[], double b[]);
void triangular_invert(int n, double a[]);
int cholesky(int n, double a[]);

void qr_insert(int n, int r, double a[], double Q[], double R[]);
void qr_delete(int n, int r, int col, double Q[], double R[]);

/*
 * code gleaned from Powell's ZQPCVX routine to determine a small
 * number that can be assumed to be an upper bound on the relative
 * precision of the computer arithmetic.
 */
double calculate_vsmall() {
    double vsmall = 1e-60;
    do {
        vsmall += vsmall;
    } while ((vsmall * .1 + 1.) <= 1. || (vsmall * .2 + 1.) <= 1.);
    return vsmall;
}

/*
 * Input parameters
 * ----------------
 * G      nxn matrix, the matrix G from above
 *        *** WILL BE DESTROYED ON EXIT ***
 * *      The user has two possibilities:
 *        a) Give G (passing factorized = 0). G must be symmetric and positive definite.
 *        b) Give R^-1 (passing factorized != 0), where R is upper triangular and G = R^T R.
 *
 *        If G is passed, then R^-1 will be computed using a generic routine.
 *        So if it is cheaper to calculate R^-1 in another way (for example if G is a band matrix),
 *        then it may be preferable to pass R^-1 directly.
 *
 * av     nx1 vector, the vector a from above
 *        *** WILL BE DESTROYED ON EXIT ***
 *
 *        On exit, contains the solution to the unconstrained problem.
 *
 * n      the dimension of G and av
 *
 * C      nxq matrix, the constraint matrix C from above (C^T = (C1 C2)^T)
 *
 * bv     qx1 vector, the constraint vector b from above (b^T = (b1 b2)^T)
 *
 * q      the number of constraints.
 *
 * meq    the number of equality constraints, 0 <= meq <= q.
 *
 * work   an array of length >= 2n + 2q + r*(r+5)/2, where r = min(n, q)
 *        for storage of intermediate values
 *
 * factorized   whether G itself or its inverted factor R^-1 is passed in the G parameter.

 * Output parameters
 * -----------------
 * xv     nx1 vector, receives the solution x to the minimisation problem
 *
 * lagr   qx1 vector, receives the Lagrange multipliers
 *
 * obj    receives the value of the objective
 *
 * iact   nx1 vector, receives in the first nact components the 1-based indices of the constraints in the active set.
 *
 * nact   receives the number of constraints in the active set.
 *
 * iter   2x1 vector:
 *        a) first component receives the number of times a constraint was added to the active set
 *           (occurs once per iteration)
 *        b) second component receives the number of times a constraint was removed from the active set
 *           (occurs zero or more times per iteration)

 * Return value
 * ------------
 * 0, solution was found
 * 1, the problem has no solution
 * 2, a matrix G was supplied that was not positive definite
 */
int qpgen2_(double *G, double *av, int n,
    double *xv, double *lagr, double *obj,
    double *C, double *bv, int q, int meq,
    int *iact, int *nact, int *iter,
    double *work, int factorized)
{
    double vsmall = calculate_vsmall();

    int* pIterFull = &iter[0];
    int* pIterPartial = &iter[1];

    int r = n <= q ? n : q;
    double *dv = work;
    double *zv = dv + n;
    double *rv = zv + n;
    double *uv = rv + r;
    double *R = uv + r;
    double *sv = R + r * (r + 1) / 2;
    double *nbv = sv + q;
    int work_length = n + n + r + r + r * (r + 1) / 2 + q + q;

    for (int i = 0; i < work_length; i++) {
        work[i] = 0.;
    }

    for (int i = 0; i < q; i++) {
        iact[i] = 0;
        lagr[i] = 0.;
    }

    // Initialisation. We want:
    // - xv and av to contain G^-1 a, the unconstrained minimum;
    // - J to contain L^-T, the inverse of the upper triangular Cholesky factor of G.

    for (int i = 0; i < n; i++) {
        xv[i] = av[i];
    }

    if (!factorized) {
        if (cholesky(n, G) != 0) { // now the upper triangle of G contains L^T
            return 2;
        }
        triangular_solve_transpose(n, G, xv); // now xv contains L^-1 a
        triangular_solve(n, G, xv);           // now xv contains L^-T L^-1 a = G^-1 a
        triangular_invert(n, G);              // now G contains L^-T
    } else {
        // G is already L^-T
        triangular_multiply_transpose(n, G, xv); // now xv contains L^-1 a
        triangular_multiply(n, G, xv);           // now xv contains L^-T L^-1 a = G^-1 a
    }

    double *J = G;

    // Set the lower triangle of J to zero.

    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < n; i++) {
            J[i + j * n] = 0.;
        }
    }

    // Calculate the objective value at the unconstrained minimum.

    *obj = -dot(n, av, xv) / 2.;

    // Store the unconstrained minimum in av for return.

    for (int i = 0; i < n; i++) {
        av[i] = xv[i];
    }

    // Calculate the norm of each column of the C matrix.
    // This will be used in our pivoting rule.

    for (int i = 0; i < q; i++) {
        nbv[i] = sqrt(dot(n, &C[i * n], &C[i * n]));
    }

    *nact = 0;
    *pIterPartial = 0;

    for (*pIterFull = 1; ; (*pIterFull)++) {

        // Calculate the slack variables C^T xv - bv and store the result in sv.

        for (int i = 0; i < q; i++) {
            double temp = dot(n, xv, &C[i * n]) - bv[i];
            sv[i] = fabs(temp) < vsmall ? 0. : temp;
        }
        // Force the slack variables to zero for constraints in the active set,
        // as a safeguard against rounding errors.
        for (int i = 0; i < *nact; i++) {
            sv[iact[i] - 1] = 0.;
        }

        // Choose a violated constraint to add to the active set.
        // We choose the constraint with the largest violation.
        // The index of the constraint to add is stored in iadd.

        int iadd = 0;
        double max_violation = 0.;
        for (int i = 0; i < q; i++) {
            if (sv[i] < -max_violation * nbv[i]) {
                iadd = i + 1;
                max_violation = -sv[i] / nbv[i];
            } else if (i < meq && sv[i] > max_violation * nbv[i]) {
                iadd = i + 1;
                max_violation = sv[i] / nbv[i];
            }
        }

        if (iadd == 0) {
            // All constraints are satisfied. We are at the optimum.

            for (int i = 0; i < *nact; i++) {
                lagr[iact[i] - 1] = uv[i];
            }
            return 0;
        }

        double slack = sv[iadd - 1];
        double reverse_step = slack > 0.;
        double u = 0;

        for (; ; (*pIterPartial)++) {
            // Set dv = J^T n, where n is the column of C corresponding to the constraint
            // that we are adding to the active set.

            for (int i = 0; i < n; i++) {
                dv[i] = dot(n, &J[i * n], &C[(iadd - 1) * n]);
            }

            // Set zv = J_2 d_2. This is the step direction for the primal variable xv.

            for (int i = 0; i < n; i++) {
                zv[i] = 0.;
            }
            for (int j = *nact; j < n; j++) {
                axpy(n, dv[j], &J[j * n], zv);
            }

            // Set rv = R^-1 d_1. This is (the negative of) the step direction for the dual variable uv.

            for (int i = 0; i < *nact; i++) {
                rv[i] = dv[i];
            }
            for (int i = *nact - 1; i >= 0; i--) {
                rv[i] /= R[(i + 1) * (i + 2) / 2 - 1];
                axpy(i, -rv[i], &R[i * (i + 1) / 2], rv);
            }

            // Find the largest step length t1 before dual feasibility is violated.
            // Store in idel the index of the constraint to remove from the active set, if we get that far.

            int t1inf = 1, idel;
            double t1;
            for (int i = 0; i < *nact; i++) {
                if (iact[i] > meq && ((!reverse_step && rv[i] > 0.) || (reverse_step && rv[i] < 0.))) {
                    double temp = uv[i] / fabs(rv[i]);
                    if (t1inf || temp < t1) {
                        t1inf = 0;
                        t1 = temp;
                        idel = i + 1;
                    }
                }
            }

            // Find the step length t2 to bring the slack variable to zero for the constraint we are adding to the active set.
            // Store in ztn the rate at which the slack variable is increased. This is used to update the objective value below.

            int t2inf = fabs(dot(n, zv, zv)) <= vsmall;
            double t2, ztn;
            if (!t2inf) {
                ztn = dot(n, zv, &C[(iadd - 1) * n]);
                t2 = fabs(slack) / ztn;
            }

            if (t1inf && t2inf) {
                // Can step infinitely far; dual problem is unbounded and primal problem is infeasible.
                return 1;
            }

            // We will take a full step if t2 <= t1.
            int full_step = !t2inf && (t1inf || t1 >= t2);
            double step_length = full_step ? t2 : t1;
            double step = reverse_step ? -step_length : step_length;

            if (!t2inf) {
                // Update primal variable
                axpy(n, step, zv, xv);

                // Update objective value
                *obj += step * ztn * (step / 2. + u);
            }

            // Update dual variable
            axpy(*nact, -step, rv, uv);
            u += step;

            if (full_step) {
                break;
            }

            // Remove constraint idel from the active set.

            qr_delete(n, *nact, idel, J, R);
            for (int i = idel; i < *nact; i++) {
                uv[i - 1] = uv[i];
                iact[i - 1] = iact[i];
            }
            uv[*nact - 1] = 0.;
            iact[*nact - 1] = 0;
            --(*nact);

            if (!t2inf) {
                // We took a step in primal space, but only took a partial step.
                // So we need to update the slack variable that we are currently bringing to zero.
                slack = dot(n, xv, &C[(iadd - 1) * n]) - bv[iadd - 1];
            }
        }

        // Add constraint iadd to the active set.

        ++(*nact);
        uv[*nact - 1] = u;
        iact[*nact - 1] = iadd;
        qr_insert(n, *nact, dv, J, R);
    }
}
