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
double hypot(double, double);
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
    } while ((vsmall * .1 + 1.) <= 1.0 || (vsmall * .2 + 1.) <= 1.0);
    return vsmall;
}

int qpgen2_(double *dmat, double *dvec, int n,
    double *sol, double *lagr, double *crval,
    double *amat, double *bvec, int q, int meq,
    int *iact, int *nact, int *iter,
    double *work, int factorized)
{
    double t1, tt;
    int it1, nvl;
    double sum;
    double temp;
    int t1inf, t2min;
    double vsmall = calculate_vsmall();

    int* pIterFull = &iter[0];
    int* pIterPartial = &iter[1];

    // calculate some constants, i.e., from which index on the different
    // quantities are stored in the work matrix

    int r = n <= q ? n : q;
    int iwzv = n;
    int iwrv = iwzv + n;
    int iwuv = iwrv + r;
    int iwrm = iwuv + r + 1;
    int iwsv = iwrm + r * (r + 1) / 2;
    int iwnbv = iwsv + q;
    double *zv = &work[iwzv];
    double *rv = &work[iwrv];
    double *uv = &work[iwuv];
    double *rm = &work[iwrm];
    double *sv = &work[iwsv];
    double *nbv = &work[iwnbv];

    for (int i = 0; i < iwnbv + q; i++) {
        work[i] = 0.;
    }

    for (int i = 0; i < q; i++) {
        iact[i] = 0;
        lagr[i] = 0.;
    }

    // get the initial solution

    for (int i = 0; i < n; i++) {
        sol[i] = dvec[i];
    }

    if (!factorized) {
        if (cholesky(n, dmat) != 0) {
            return 2;
        }
        triangular_solve_transpose(n, dmat, sol);
        triangular_solve(n, dmat, sol);
        triangular_invert(n, dmat);
    } else {
        // Matrix D is already factorized, so we have to multiply d first with
        // R^-T and then with R^-1.  R^-1 is stored in the upper half of the
        // array dmat.

        for (int j = n - 1; j >= 0; j--) {
            sol[j] *= dmat[j + j * n];
            for (int i = 0; i < j; i++) {
                sol[j] += dmat[i + j * n] * sol[i];
            }
        }

        for (int j = 0; j < n; j++) {
            sol[j] *= dmat[j + j * n];
            for (int i = j + 1; i < n; i++) {
                sol[j] += dmat[j + i * n] * sol[i];
            }
        }
    }

    // set lower triangular of dmat to zero

    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < n; i++) {
            dmat[i + j * n] = 0.;
        }
    }

    // calculate value of the criterion at unconstrained minima

    *crval = -dot(n, dvec, sol) / 2.;

    // now we can return the unconstrained minimum in dvec

    for (int i = 0; i < n; i++) {
        dvec[i] = sol[i];
    }

    // calculate the norm of each column of the A matrix

    for (int i = 0; i < q; i++) {
        nbv[i] = sqrt(dot(n, &amat[i * n], &amat[i * n]));
    }

    *nact = 0;
    *pIterFull = 0;
    *pIterPartial = 0;

DONE_FULL_STEP:

    // start a new iteration

    ++*pIterFull;

    // calculate all constraints and check which are still violated
    // for the equality constraints we have to check whether the normal
    // vector has to be negated (as well as bvec in that case)

    for (int i = 0; i < q; i++) {
        sum = -bvec[i];
        for (int j = 0; j < n; j++) {
            sum += amat[j + i * n] * sol[j];
        }
        if (fabs(sum) < vsmall) {
            sum = 0.;
        }
        if (i >= meq) {
            sv[i] = sum;
        } else {
            sv[i] = -fabs(sum);
            if (sum > 0.) {
                for (int j = 0; j < n; j++) {
                    amat[j + i * n] = -amat[j + i * n];
                }
                bvec[i] = -bvec[i];
            }
        }
    }

    // as safeguard against rounding errors set already active constraints
    // explicitly to zero

    for (int i = 0; i < *nact; i++) {
        sv[iact[i] - 1] = 0.;
    }

    // we weight each violation by the number of non-zero elements in the
    // corresponding row of A. then we choose the violated constraint which
    // has maximal absolute value, i.e., the minimum.

    nvl = 0;
    temp = 0.;
    for (int i = 0; i < q; i++) {
        if (sv[i] < temp * nbv[i]) {
            nvl = i + 1;
            temp = sv[i] / nbv[i];
        }
    }
    if (nvl == 0) {
        for (int i = 0; i < *nact; i++) {
            lagr[iact[i] - 1] = uv[i];
        }
        return 0;
    }

DONE_PARTIAL_STEP:

    // calculate d=J^Tn^+ where n^+ is the normal vector of the violated
    // constraint. J is stored in dmat in this implementation!!
    // if we drop a constraint, we have to jump back here.

    for (int i = 0; i < n; i++) {
        sum = 0.;
        for (int j = 0; j < n; j++) {
            sum += dmat[j + i * n] * amat[j + (nvl - 1) * n];
        }
        work[i] = sum;
    }

    // Now calculate z = J_2 d_2

    for (int i = 0; i < n; i++) {
        zv[i] = 0.;
    }
    for (int j = *nact; j < n; j++) {
        for (int i = 0; i < n; i++) {
            zv[i] += dmat[i + j * n] * work[j];
        }
    }

    // and r = R^{-1} d_1, check also if r has positive elements (among the
    // entries corresponding to inequalities constraints).

    for (int i = *nact - 1; i >= 0; i--) {
        sum = work[i];
        int rm_offset = (i + 1) * (i + 2) / 2 - 1;
        int k = 0;
        for (int j = i + 1; j <= *nact; j++) {
            sum -= rm[rm_offset + i + 1 + k] * rv[j];
            k += j + 1;
        }
        rv[i] = sum / rm[rm_offset];
    }

    // if r has positive elements, find the partial step length t1, which is
    // the maximum step in dual space without violating dual feasibility.
    // it1  stores in which component t1, the min of u/r, occurs.

    t1inf = 1;
    for (int i = 0; i < *nact; i++) {
        if (iact[i] > meq && rv[i] > 0.) {
            double current = uv[i] / rv[i];
            if (t1inf) {
                t1inf = 0;
                t1 = current;
                it1 = i + 1;
            } else if (current < t1) {
                it1 = i + 1;
            }
        }
    }

    // test if the z vector is equal to zero

    sum = 0.;
    for (int i = 0; i < n; i++) {
        sum += zv[i] * zv[i];
    }
    if (fabs(sum) <= vsmall) {
        // No step in primal space such that the new constraint becomes
        // feasible. Take step in dual space and drop a constant.

        if (t1inf) {
            // No step in dual space possible either, problem is not solvable
            return 1;
        }

        // we take a partial step in dual space and drop constraint it1,
        // that is, we drop the it1-th active constraint.
        // then we continue at step 2(a) (marked by label 55)

        for (int i = 0; i < *nact; i++) {
            uv[i] -= t1 * rv[i];
        }
        uv[*nact] += t1;
    } else {
        // compute full step length t2, minimum step in primal space such that */
        // the constraint becomes feasible. */
        // keep sum (which is z^Tn^+) to update crval below! */

        sum = 0.;
        for (int i = 0; i < n; i++) {
            sum += zv[i] * amat[i + (nvl - 1) * n];
        }
        tt = -sv[nvl - 1] / sum;
        t2min = 1;
        if (!t1inf && t1 < tt) {
            tt = t1;
            t2min = 0;
        }

        // take step in primal and dual space

        for (int i = 0; i < n; i++) {
            sol[i] += tt * zv[i];
        }
        *crval += tt * sum * (tt / 2. + uv[*nact]);
        for (int i = 0; i < *nact; i++) {
            uv[i] -= tt * rv[i];
        }
        uv[*nact] += tt;

        // if it was a full step, then we check wheter further constraints are
        // violated otherwise we can drop the current constraint and iterate once
        // more
        if (t2min) {
            // we took a full step. Thus add constraint nvl to the list of active
            // constraints and update J and R

            ++(*nact);
            iact[*nact - 1] = nvl;
            qr_insert(n, *nact, work, dmat, rm);

            goto DONE_FULL_STEP;
        } else {
            // we took a partial step in dual space. Thus drop constraint it1,
            // that is, we drop the it1-th active constraint.
            // then we continue at step 2(a) (marked by label 55)
            // but since the fit changed, we have to recalculate now "how much"
            // the fit violates the chosen constraint now.

            sum = -bvec[nvl-1];
            for (int j = 0; j < n; j++) {
                sum += sol[j] * amat[j + (nvl-1) * n];
            }
            if (nvl > meq) {
                sv[nvl - 1] = sum;
            } else {
                sv[nvl - 1] = -fabs(sum);
                if (sum > 0.) {
                    for (int j = 0; j < n; j++) {
                        amat[j + (nvl - 1) * n] = -amat[j + (nvl - 1) * n];
                    }
                    bvec[nvl - 1] = -bvec[nvl - 1];
                }
            }
        }
    }

    // Drop constraint it1
    qr_delete(n, *nact, it1, dmat, rm);
    for (int i = it1; i < *nact; i++) {
        uv[i - 1] = uv[i];
        iact[i - 1] = iact[i];
    }
    uv[*nact-1] = uv[*nact];
    uv[*nact] = 0.;
    iact[*nact-1] = 0;
    --(*nact);
    ++(*pIterPartial);

    goto DONE_PARTIAL_STEP;
}
