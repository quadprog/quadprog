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

void triangular_solve(int n, double a[], double b[]);
void triangular_solve_transpose(int n, double a[], double b[]);
void triangular_invert(int n, double a[]);
int cholesky(int n, double a[]);

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

    int* pIterMain = &iter[0];
    int* pIterDeleted = &iter[1];

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

    // store the initial dvec to calculate below the unconstrained minima of
    // the critical value.

    for (int i = 0; i < n; i++) {
        work[i] = dvec[i];
    }
    for (int i = n; i < iwnbv + q; i++) {
        work[i] = 0.;
    }
    for (int i = 0; i < q; i++) {
        iact[i] = 0;
        lagr[i] = 0.;
    }

    // get the initial solution

    if (!factorized) {
        if (cholesky(n, dmat) != 0) {
            return 2;
        }
        triangular_solve_transpose(n, dmat, dvec);
        triangular_solve(n, dmat, dvec);
        triangular_invert(n, dmat);
    } else {
        // Matrix D is already factorized, so we have to multiply d first with
        // R^-T and then with R^-1.  R^-1 is stored in the upper half of the
        // array dmat.

        for (int j = 0; j < n; j++) {
            sol[j] = 0.;
            for (int i = 0; i <= j; i++) {
                sol[j] += dmat[i + j * n] * dvec[i];
            }
        }

        for (int j = 0; j < n; j++) {
            dvec[j] = 0.;
            for (int i = j; i < n; i++) {
                dvec[j] += dmat[j + i * n] * sol[i];
            }
        }
    }

    // set lower triangular of dmat to zero, store dvec in sol and
    // calculate value of the criterion at unconstrained minima

    *crval = 0.;
    for (int j = 0; j < n; j++) {
        sol[j] = dvec[j];
        *crval += work[j] * sol[j];
        work[j] = 0.;
        for (int i = j + 1; i < n; i++) {
            dmat[i + j * n] = 0.;
        }
    }
    *crval = -(*crval) / 2.;

    // calculate the norm of each column of the A matrix

    for (int i = 0; i < q; i++) {
        sum = 0.;
        for (int j = 0; j < n; j++) {
            sum += amat[j + i * n] * amat[j + i * n];
        }
        nbv[i] = sqrt(sum);
    }

    *nact = 0;
    *pIterMain = 0;
    *pIterDeleted = 0;

DONE_FULL_STEP:

    // start a new iteration

    ++*pIterMain;

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

            // to update R we have to put the first nact-1 components of the d vector
            // into column (nact) of R

            for (int i = 0; i < *nact - 1; i++) {
                rm[(*nact - 1) * *nact / 2 + i] = work[i];
            }

            // if now nact=n, then we just have to add the last element to the new
            // row of R.
            // Otherwise we use Givens transformations to turn the vector d(nact:n)
            // into a multiple of the first unit vector. That multiple goes into the
            // last element of the new row of R and J is accordingly updated by the
            // Givens transformations.

            for (int i = n - 1; i >= *nact; i--) {
                // we have to find the Givens rotation which will reduce the element
                // (l1) of d to zero.
                // if it is already zero we don't have to do anything, except of
                // decreasing l1

                if (work[i] != 0. && work[i - 1] == 0.) {
                    work[i - 1] = work[i];
                    for (int j = 0; j < n; j++) {
                        temp = dmat[j + (i - 1) * n];
                        dmat[j + (i - 1) * n] = dmat[j + i * n];
                        dmat[j + i * n] = temp;
                    }
                } else if (work[i] != 0.) {
                    double h = hypot(work[i - 1], work[i]);

                    if (work[i - 1] < 0.0) {
                        h = -h;
                    }

                    double gc = work[i - 1] / h;
                    double gs = work[i] / h;
                    double nu = work[i] / (work[i - 1] + h);

                    // The Givens rotation is done with the matrix (gc gs, gs -gc).
                    // If gc is one, then element (i) of d is zero compared with element
                    // (l1-1). Hence we don't have to do anything.
                    // If gc is zero, then we just have to switch column (i) and column (i-1)
                    // of J. Since we only switch columns in J, we have to be careful how we
                    // update d depending on the sign of gs.
                    // Otherwise we have to apply the Givens rotation to these columns.
                    // The i-1 element of d has to be updated to temp.

                    work[i - 1] = h;
                    for (int j = 0; j < n; j++) {
                        temp = gc * dmat[j + (i - 1) * n] + gs * dmat[j + i * n];
                        dmat[j + i * n] = nu * (dmat[j + (i - 1) * n] + temp) - dmat[j + i * n];
                        dmat[j + (i - 1) * n] = temp;
                    }
                }
            }

            // So store d(nact) in R(nact,nact)
            rm[(*nact) * (*nact + 1) / 2 - 1] = work[*nact - 1];

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

    // if it1 = nact it is only necessary to update the vector u and nact
    for (; it1 < *nact; it1++) {
        // After updating one row of R (column of J) we will also come back here

        // we have to find the Givens rotation which will reduce the element
        // (it1+1,it1+1) of R to zero.
        // if it is already zero we don't have to do anything except of updating
        // u, iact, and shifting column (it1+1) of R to column (it1)
        // l  will point to element (1,it1+1) of R
        // l1 will point to element (it1+1,it1+1) of R

        int l = it1 * (it1 + 1) / 2;
        int l1 = l + it1;
        if (rm[l1] != 0. && rm[l1 - 1] == 0.) {
            for (int i = it1 + 1; i <= *nact; i++) {
                temp = rm[l1 - 1];
                rm[l1 - 1] = rm[l1];
                rm[l1] = temp;
                l1 += i;
            }
            for (int i = 0; i < n; i++) {
                temp = dmat[i + (it1-1) * n];
                dmat[i + (it1 - 1) * n] = dmat[i + (it1) * n];
                dmat[i + (it1) * n] = temp;
            }
        } else if (rm[l1] != 0.) {
            double h = hypot(rm[l1 - 1], rm[l1]);

            if (rm[l1 - 1] < 0.0) {
                h = -h;
            }

            double gc = rm[l1 - 1] / h;
            double gs = rm[l1] / h;
            double nu = rm[l1] / (rm[l1 - 1] + h);

            // The Givens rotatin is done with the matrix (gc gs, gs -gc).
            // If gc is one, then element (it1+1,it1+1) of R is zero compared with
            // element (it1,it1+1). Hence we don't have to do anything.
            // if gc is zero, then we just have to switch row (it1) and row (it1+1)
            // of R and column (it1) and column (it1+1) of J. Since we swithc rows in
            // R and columns in J, we can ignore the sign of gs.
            // Otherwise we have to apply the Givens rotation to these rows/columns.

            for (int i = it1 + 1; i <= *nact; i++) {
                temp = gc * rm[l1 - 1] + gs * rm[l1];
                rm[l1] = nu * (rm[l1 - 1] + temp) - rm[l1];
                rm[l1 - 1] = temp;
                l1 += i;
            }
            for (int i = 0; i < n; i++) {
                temp = gc * dmat[i + (it1-1) * n] + gs * dmat[i + (it1) * n];
                dmat[i + (it1) * n] = nu * (dmat[i + (it1-1) * n] + temp) - dmat[i + (it1) * n];
                dmat[i + (it1-1) * n] = temp;
            }
        }

        // shift column (it1+1) of R to column (it1) (that is, the first it1
        // elements). The posit1on of element (1,it1+1) of R was calculated above
        // and stored in l.

        for (int i = 0; i < it1; i++) {
            rm[l + i - it1] = rm[l + i];
        }

        // update vector u and iact as necessary
        // Continue with updating the matrices J and R

        uv[it1-1] = uv[it1];
        iact[it1-1] = iact[it1];
    }

    uv[*nact-1] = uv[*nact];
    uv[*nact] = 0.;
    iact[*nact-1] = 0;
    --(*nact);
    ++(*pIterDeleted);

    goto DONE_PARTIAL_STEP;
}
