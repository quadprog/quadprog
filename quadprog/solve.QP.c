/*  Copyright (C) 1995-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com> */

/*  This program is free software; you can redistribute it and/or modify */
/*  it under the terms of the GNU General Public License as published by */
/*  the Free Software Foundation; either version 2 of the License, or */
/*  (at your option) any later version. */

/*  This program is distributed in the hope that it will be useful, */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*  GNU General Public License for more details. */

/*  You should have received a copy of the GNU General Public License */
/*  along with this program; if not, write to the Free Software */
/*  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, */
/*  USA. */

/*  this routine uses the Goldfarb/Idnani algorithm to solve the */
/*  following minimization problem: */

/*        minimize  -d^T x + 1/2 *  x^T D x */
/*        where   A1^T x  = b1 */
/*                A2^T x >= b2 */

/*  the matrix D is assumed to be positive definite.  Especially, */
/*  w.l.o.g. D is assumed to be symmetric. */

/*  Input parameter: */
/*  dmat   nxn matrix, the matrix D from above (dp) */
/*         *** WILL BE DESTROYED ON EXIT *** */
/*         The user has two possibilities: */
/*         a) Give D (factorized=0), in this case we use routines from LINPACK */
/*            to decompose D. */
/*         b) To get the algorithm started we need R^-1, where D=R^TR. */
/*            So if it is cheaper to calculate R^-1 in another way (D may */
/*            be a band matrix) then with the general routine, the user */
/*            may pass R^{-1}.  Indicated by factorized not equal to zero. */
/*  dvec   nx1 vector, the vector d from above (dp) */
/*         *** WILL BE DESTROYED ON EXIT *** */
/*         contains on exit the solution to the initial, i.e., */
/*         unconstrained problem */
/*  n      the dimension of dmat and dvec (int) */
/*  amat   nxq matrix, the matrix A from above (dp) [ A=(A1 A2)^T ] */
/*         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE */
/*             CHANGED SIGNES ON EXIT *** */
/*  bvec   qx1 vector, the vector of constants b in the constraints (dp) */
/*         [ b = (b1^T b2^T)^T ] */
/*         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE */
/*             CHANGED SIGNES ON EXIT *** */
/*  q      integer, the number of constraints. */
/*  meq    integer, the number of equality constraints, 0 <= meq <= q. */
/*  factorized   integer, code for the status of the matrix D: */
/*            factorized =  0, we have to decompose D */
/*            factorized != 0, D is already decomposed into D=R^TR and we were */
/*                       given R^{-1}. */

/*  Output parameter: */
/*  sol   nx1 the final solution (x in the notation above) */
/*  lagr  qx1 the final Lagrange multipliers */
/*  crval scalar, the value of the criterion at the minimum */
/*  iact  qx1 vector, the constraints which are active in the final */
/*        fit (int) */
/*  nact  scalar, the number of constraints active in the final fit (int) */
/*  iter  2x1 vector, first component gives the number of "main" */
/*        iterations, the second one says how many constraints were */
/*        deleted after they became active */

/*  Return value: */
/*        0, no problems */
/*        1, the minimization problem has no solution */
/*        2, problems with decomposing D, in this case sol */
/*           contains garbage!! */

/*  Working space: */
/*  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1 */
/*        where r=min(n,q) */

double fabs(double);
double sqrt(double);

double dot(int n, double x[], double y[]);
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

int qpgen2_(double *dmat, double *dvec, int n,
    double *sol, double *lagr, double *crval,
    double *amat, double *bvec, int q, int meq,
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
    double *rm = uv + r + 1;
    double *sv = rm + r * (r + 1) / 2;
    double *nbv = sv + q;
    int work_length = n + n + r + (r + 1) + r * (r + 1) / 2 + q + q;

    for (int i = 0; i < work_length; i++) {
        work[i] = 0.;
    }

    for (int i = 0; i < q; i++) {
        iact[i] = 0;
        lagr[i] = 0.;
    }

    // Initialisation. We want:
    // - sol and dvec to contain G^-1 a, the unconstrained minimum;
    // - dmat to contain L^-T, the inverse of the upper triangular Cholesky factor of G.

    for (int i = 0; i < n; i++) {
        sol[i] = dvec[i];
    }

    if (!factorized) {
        if (cholesky(n, dmat) != 0) { // now the upper triangle of dmat contains L^T
            return 2;
        }
        triangular_solve_transpose(n, dmat, sol); // now sol contains L^-1 a
        triangular_solve(n, dmat, sol);           // now sol contains L^-T L^-1 a = G^-1 a
        triangular_invert(n, dmat);               // now dmat contains L^-T
    } else {
        // dmat is already L^-T

        // multiply sol by L^-T
        for (int j = n - 1; j >= 0; j--) {
            sol[j] *= dmat[j + j * n];
            for (int i = 0; i < j; i++) {
                sol[j] += dmat[i + j * n] * sol[i];
            }
        }

        // multiply sol by L^-1
        for (int j = 0; j < n; j++) {
            sol[j] *= dmat[j + j * n];
            for (int i = j + 1; i < n; i++) {
                sol[j] += dmat[j + i * n] * sol[i];
            }
        }
    }

    // Set the lower triangle of dmat to zero.

    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < n; i++) {
            dmat[i + j * n] = 0.;
        }
    }

    // Calculate the objective value at the unconstrained minimum.

    *crval = -dot(n, dvec, sol) / 2.;

    // Store the unconstrained minimum in dvec for return.

    for (int i = 0; i < n; i++) {
        dvec[i] = sol[i];
    }

    // Calculate the norm of each column of the A matrix.
    // This will be used in our pivoting rule.

    for (int i = 0; i < q; i++) {
        nbv[i] = sqrt(dot(n, &amat[i * n], &amat[i * n]));
    }

    *nact = 0;
    *pIterPartial = 0;

    for (*pIterFull = 1; ; (*pIterFull)++) {

        // Calculate the slack variables A^T x - b and store the result in sv.

        for (int i = 0; i < q; i++) {
            double temp = dot(n, sol, &amat[i * n]) - bvec[i];
            if (fabs(temp) < vsmall) {
                temp = 0.;
            }
            if (i >= meq) {
                sv[i] = temp;
            } else {
                sv[i] = -fabs(temp);
                if (temp > 0.) {
                    for (int j = 0; j < n; j++) {
                        amat[j + i * n] = -amat[j + i * n];
                    }
                    bvec[i] = -bvec[i];
                }
            }
        }
        // Force the slack variables to zero for constraints in the active set,
        // as a safeguard against rounding errors.
        for (int i = 0; i < *nact; i++) {
            sv[iact[i] - 1] = 0.;
        }

        // Choose a violated constraint to add to the active set.
        // We choose the constraint with the largest violation.
        // The index of the constraint to add is stored in nvl.

        int nvl = 0;
        double max_violation = 0.;
        for (int i = 0; i < q; i++) {
            if (sv[i] < max_violation * nbv[i]) {
                nvl = i + 1;
                max_violation = sv[i] / nbv[i];
            }
        }

        if (nvl == 0) {
            // All constraints are satisfied. We are at the optimum.

            for (int i = 0; i < *nact; i++) {
                lagr[iact[i] - 1] = uv[i];
            }
            return 0;
        }

        double slack = sv[nvl - 1];

        for (; ; (*pIterPartial)++) {
            // Set dv = J^T n, where n is the column of A corresponding to the constraint
            // that we are adding to the active set.

            for (int i = 0; i < n; i++) {
                dv[i] = dot(n, &dmat[i * n], &amat[(nvl - 1) * n]);
            }

            // Set zv = J_2 d_2. This is the step direction for the primal variable sol.

            for (int i = 0; i < n; i++) {
                zv[i] = 0.;
            }
            for (int j = *nact; j < n; j++) {
                axpy(n, dv[j], &dmat[j * n], zv);
            }

            // Set rv = R^-1 d_1. This is (the negative of) the step direction for the dual variable uv.

            for (int i = *nact - 1; i >= 0; i--) {
                double temp = dv[i];
                int k = 0;
                for (int j = i + 1; j < *nact; j++) {
                    temp -= rm[(i + 2) * (i + 3) / 2 - 2 + k] * rv[j];
                    k += j + 1;
                }
                rv[i] = temp / rm[(i + 1) * (i + 2) / 2 - 1];
            }

            // Find the largest step length t1 before dual feasibility is violated.
            // Store in it1 the index of the constraint to remove from the active set, if we get that far.

            int t1inf = 1, it1;
            double t1;
            for (int i = 0; i < *nact; i++) {
                if (iact[i] > meq && rv[i] > 0.) {
                    double temp = uv[i] / rv[i];
                    if (t1inf) {
                        t1inf = 0;
                        t1 = temp;
                        it1 = i + 1;
                    } else if (temp < t1) {
                        it1 = i + 1;
                    }
                }
            }

            // Find the step length t2 to bring the slack variable to zero for the constraint we are adding to the active set.
            // Store in ztn the rate at which the slack variable is increased. This is used to update the objective value below.

            int t2inf = fabs(dot(n, zv, zv)) <= vsmall;
            double t2, ztn;
            if (!t2inf) {
                ztn = dot(n, zv, &amat[(nvl - 1) * n]);
                t2 = -slack / ztn;
            }

            if (t1inf && t2inf) {
                // Can step infinitely far; dual problem is unbounded and primal problem is infeasible.
                return 1;
            }

            // We will take a full step if t2 <= t1.
            int t2min = !t2inf && (t1inf || t1 >= t2);
            double tt = t2min ? t2 : t1;

            if (!t2inf) {
                // Update primal variable
                axpy(n, tt, zv, sol);

                // Update objective value
                *crval += tt * ztn * (tt / 2. + uv[*nact]);
            }

            // Update dual variable
            axpy(*nact, -tt, rv, uv);
            uv[*nact] += tt;

            if (t2min) {
                break;
            }

            // Remove constraint it1 from the active set.

            qr_delete(n, *nact, it1, dmat, rm);
            for (int i = it1; i < *nact; i++) {
                uv[i - 1] = uv[i];
                iact[i - 1] = iact[i];
            }
            uv[*nact - 1] = uv[*nact];
            uv[*nact] = 0.;
            iact[*nact-1] = 0;
            --(*nact);

            if (!t2inf) {
                // We took a step in primal space, but only took a partial step.
                // So we need to update the slack variable that we are currently bringing to zero.

                double temp = dot(n, sol, &amat[(nvl - 1) * n]) - bvec[nvl - 1];
                if (nvl > meq) {
                    slack = temp;
                } else {
                    slack = -fabs(temp);
                    if (temp > 0.) {
                        for (int j = 0; j < n; j++) {
                            amat[j + (nvl - 1) * n] = -amat[j + (nvl - 1) * n];
                        }
                        bvec[nvl - 1] = -bvec[nvl - 1];
                    }
                }
            }
        }

        // Add constraint nvl to the active set.

        ++(*nact);
        iact[*nact - 1] = nvl;
        qr_insert(n, *nact, dv, dmat, rm);
    }
}
