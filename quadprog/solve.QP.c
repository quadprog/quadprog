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
/*         a) Give D (ierr=0), in this case we use routines from LINPACK */
/*            to decompose D. */
/*         b) To get the algorithm started we need R^-1, where D=R^TR. */
/*            So if it is cheaper to calculate R^-1 in another way (D may */
/*            be a band matrix) then with the general routine, the user */
/*            may pass R^{-1}.  Indicated by ierr not equal to zero. */
/*  dvec   nx1 vector, the vector d from above (dp) */
/*         *** WILL BE DESTROYED ON EXIT *** */
/*         contains on exit the solution to the initial, i.e., */
/*         unconstrained problem */
/*  fddmat scalar, the leading dimension of the matrix dmat */
/*  n      the dimension of dmat and dvec (int) */
/*  amat   nxq matrix, the matrix A from above (dp) [ A=(A1 A2)^T ] */
/*         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE */
/*             CHANGED SIGNES ON EXIT *** */
/*  bvec   qx1 vector, the vector of constants b in the constraints (dp) */
/*         [ b = (b1^T b2^T)^T ] */
/*         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE */
/*             CHANGED SIGNES ON EXIT *** */
/*  fdamat the first dimension of amat as declared in the calling program. */
/*         fdamat >= n !! */
/*  q      integer, the number of constraints. */
/*  meq    integer, the number of equality constraints, 0 <= meq <= q. */
/*  ierr   integer, code for the status of the matrix D: */
/*            ierr =  0, we have to decompose D */
/*            ierr != 0, D is already decomposed into D=R^TR and we were */
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
/*  ierr  integer, error code on exit, if */
/*           ierr = 0, no problems */
/*           ierr = 1, the minimization problem has no solution */
/*           ierr = 2, problems with decomposing D, in this case sol */
/*                     contains garbage!! */

/*  Working space: */
/*  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1 */
/*        where r=min(n,q) */

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

double d_sign(double *a, double *b) {
	double x = (*a >= 0 ? *a : -*a);
	return *b >= 0 ? x : -x;
}

/* Subroutine */ int qpgen2_(double *dmat, double *dvec, int *
	fddmat, int *n, double *sol, double *lagr, double *
	crval, double *amat, double *bvec, int *fdamat, int *
	q, int *meq, int *iact, int *nact, int *iter,
	double *work, int *ierr)
{
    /* System generated locals */
    int dmat_dim1, amat_dim1;
    double d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    int l, r__, l1;
    double t1, gc, gs, nu, tt;
    int it1, nvl;
    double sum;
    int info;
    double tmpa, tmpb, temp;
    int iwrm, iwrv, iwsv, iwuv, iwzv;
    int t1inf, t2min;
    int iwnbv;
    double vsmall;

    dmat_dim1 = *fddmat;
    amat_dim1 = *fdamat;

    int* pIterMain = &iter[0];
	int* pIterDeleted = &iter[1];

    /* Function Body */
    r__ = min(*n,*q);
    l = (*n << 1) + r__ * (r__ + 5) / 2 + (*q << 1) + 1;

/*     code gleaned from Powell's ZQPCVX routine to determine a small */
/*     number  that can be assumed to be an upper bound on the relative */
/*     precision of the computer arithmetic. */

    vsmall = 1e-60;
L1:
    vsmall += vsmall;
    tmpa = vsmall * .1 + 1.;
    tmpb = vsmall * .2 + 1.;
    if (tmpa <= 1.) {
	goto L1;
    }
    if (tmpb <= 1.) {
	goto L1;
    }

/* store the initial dvec to calculate below the unconstrained minima of */
/* the critical value. */

    for (int i = 0; i < *n; i++) {
	work[i] = dvec[i];
/* L10: */
    }
    for (int i = *n; i < l; i++) {
	work[i] = 0.;
/* L11: */
    }
    for (int i = 0; i < *q; i++) {
	iact[i] = 0;
	lagr[i] = 0.;
/* L12: */
    }

/* get the initial solution */

    if (*ierr == 0) {
	info = cholesky(*n, dmat);
	if (info != 0) {
	    *ierr = 2;
	    goto L999;
	}
	triangular_solve_transpose(*n, dmat, dvec);
	triangular_solve(*n, dmat, dvec);
	triangular_invert(*n, dmat);
    } else {

/* Matrix D is already factorized, so we have to multiply d first with */
/* R^-T and then with R^-1.  R^-1 is stored in the upper half of the */
/* array dmat. */

	for (int j = 0; j < *n; j++) {
	    sol[j] = 0.;
	    for (int i = 0; i <= j; i++) {
		sol[j] += dmat[i + j * dmat_dim1] * dvec[i];
/* L21: */
	    }
/* L20: */
	}
	for (int j = 0; j < *n; j++) {
	    dvec[j] = 0.;
	    for (int i = j; i < *n; i++) {
		dvec[j] += dmat[j + i * dmat_dim1] * sol[i];
/* L23: */
	    }
/* L22: */
	}
    }

/* set lower triangular of dmat to zero, store dvec in sol and */
/* calculate value of the criterion at unconstrained minima */

    *crval = 0.;
    for (int j = 0; j < *n; j++) {
	sol[j] = dvec[j];
	*crval += work[j] * sol[j];
	work[j] = 0.;
	for (int i = j + 1; i < *n; i++) {
	    dmat[i + j * dmat_dim1] = 0.;
/* L32: */
	}
/* L30: */
    }
    *crval = -(*crval) / 2.;
    *ierr = 0;

/* calculate some constants, i.e., from which index on the different */
/* quantities are stored in the work matrix */

    iwzv = *n;
    iwrv = iwzv + *n;
    iwuv = iwrv + r__;
    iwrm = iwuv + r__ + 1;
    iwsv = iwrm + r__ * (r__ + 1) / 2;
    iwnbv = iwsv + *q;

/* calculate the norm of each column of the A matrix */

    for (int i = 0; i < *q; i++) {
	sum = 0.;
	for (int j = 0; j < *n; j++) {
	    sum += amat[j + i * amat_dim1] * amat[j + i * amat_dim1];
/* L52: */
	}
	work[iwnbv + i] = sqrt(sum);
/* L51: */
    }
    *nact = 0;
    *pIterMain = 0;
    *pIterDeleted = 0;
L50:

/* start a new iteration */

    ++*pIterMain;

/* calculate all constraints and check which are still violated */
/* for the equality constraints we have to check whether the normal */
/* vector has to be negated (as well as bvec in that case) */

    for (int i = 0; i < *q; i++) {
	sum = -bvec[i];
	for (int j = 0; j < *n; j++) {
	    sum += amat[j + i * amat_dim1] * sol[j];
/* L61: */
	}
	if (abs(sum) < vsmall) {
	    sum = 0.;
	}
	if (i >= *meq) {
	    work[iwsv + i] = sum;
	} else {
	    work[iwsv + i] = -abs(sum);
	    if (sum > 0.) {
		for (int j = 0; j < *n; j++) {
		    amat[j + i * amat_dim1] = -amat[j + i * amat_dim1];
/* L62: */
		}
		bvec[i] = -bvec[i];
	    }
	}
/* L60: */
    }

/* as safeguard against rounding errors set already active constraints */
/* explicitly to zero */

    for (int i = 0; i < *nact; i++) {
	work[iwsv + iact[i] - 1] = 0.;
/* L70: */
    }

/* we weight each violation by the number of non-zero elements in the */
/* corresponding row of A. then we choose the violated constraint which */
/* has maximal absolute value, i.e., the minimum. */


    nvl = 0;
    temp = 0.;
    for (int i = 0; i < *q; i++) {
	if (work[iwsv + i] < temp * work[iwnbv + i]) {
	    nvl = i + 1;
	    temp = work[iwsv + i] / work[iwnbv + i];
	}
/* L71: */
    }
/* L72: */
    if (nvl == 0) {
	for (int i = 0; i < *nact; i++) {
	    lagr[iact[i] - 1] = work[iwuv + i];
/* L73: */
	}
	goto L999;
    }

/* calculate d=J^Tn^+ where n^+ is the normal vector of the violated */
/* constraint. J is stored in dmat in this implementation!! */
/* if we drop a constraint, we have to jump back here. */

L55:
    for (int i = 0; i < *n; i++) {
	sum = 0.;
	for (int j = 0; j < *n; j++) {
	    sum += dmat[j + i * dmat_dim1] * amat[j + (nvl - 1) * amat_dim1];
/* L81: */
	}
	work[i] = sum;
/* L80: */
    }

/* Now calculate z = J_2 d_2 */

    l1 = iwzv;
    for (int i = 0; i < *n; i++) {
	work[iwzv + i] = 0.;
/* L90: */
    }
    for (int j = *nact; j < *n; j++) {
	for (int i = 0; i < *n; i++) {
	    work[iwzv + i] += dmat[i + j * dmat_dim1] * work[j];
/* L93: */
	}
/* L92: */
    }

/* and r = R^{-1} d_1, check also if r has positive elements (among the */
/* entries corresponding to inequalities constraints). */

    t1inf = 1;
    for (int i = *nact - 1; i >= 0; i--) {
	sum = work[i];
	l = iwrm + (i + 1) * (i + 4) / 2 - 1;
	l1 = l - i - 1;
	for (int j = i + 1; j <= *nact; j++) {
	    sum -= work[l] * work[iwrv + j];
	    l += j + 1;
/* L96: */
	}
	sum /= work[l1];
	work[iwrv + i] = sum;
	if (iact[i] <= *meq) {
	    goto L95;
	}
	if (sum <= 0.) {
	    goto L95;
	}
/* L7: */
	t1inf = 0;
	it1 = i + 1;
L95:
	;
    }

/* if r has positive elements, find the partial step length t1, which is */
/* the maximum step in dual space without violating dual feasibility. */
/* it1  stores in which component t1, the min of u/r, occurs. */

    if (! t1inf) {
	t1 = work[iwuv + it1-1] / work[iwrv + it1-1];
	for (int i = 0; i < *nact; i++) {
	    if (iact[i] <= *meq) {
		goto L100;
	    }
	    if (work[iwrv + i] <= 0.) {
		goto L100;
	    }
	    temp = work[iwuv + i] / work[iwrv + i];
	    if (temp < t1) {
		t1 = temp;
		it1 = i + 1;
	    }
L100:
	    ;
	}
    }

/* test if the z vector is equal to zero */

    sum = 0.;
    for (int i = 0; i < *n; i++) {
	sum += work[iwzv + i] * work[iwzv + i];
/* L110: */
    }
    if (abs(sum) <= vsmall) {

/* No step in primal space such that the new constraint becomes */
/* feasible. Take step in dual space and drop a constant. */

	if (t1inf) {

/* No step in dual space possible either, problem is not solvable */

	    *ierr = 1;
	    goto L999;
	} else {

/* we take a partial step in dual space and drop constraint it1, */
/* that is, we drop the it1-th active constraint. */
/* then we continue at step 2(a) (marked by label 55) */

	    for (int i = 0; i < *nact; i++) {
		work[iwuv + i] -= t1 * work[iwrv + i];
/* L111: */
	    }
	    work[iwuv + *nact] += t1;
	    goto L700;
	}
    } else {

/* compute full step length t2, minimum step in primal space such that */
/* the constraint becomes feasible. */
/* keep sum (which is z^Tn^+) to update crval below! */

	sum = 0.;
	for (int i = 0; i < *n; i++) {
	    sum += work[iwzv + i] * amat[i + (nvl-1) * amat_dim1];
/* L120: */
	}
	tt = -work[iwsv + nvl-1] / sum;
	t2min = 1;
	if (! t1inf) {
	    if (t1 < tt) {
		tt = t1;
		t2min = 0;
	    }
	}

/* take step in primal and dual space */

	for (int i = 0; i < *n; i++) {
	    sol[i] += tt * work[iwzv + i];
/* L130: */
	}
	*crval += tt * sum * (tt / 2. + work[iwuv + *nact]);
	for (int i = 0; i < *nact; i++) {
	    work[iwuv + i] -= tt * work[iwrv + i];
/* L131: */
	}
	work[iwuv + *nact] += tt;

/* if it was a full step, then we check wheter further constraints are */
/* violated otherwise we can drop the current constraint and iterate once */
/* more */
	if (t2min) {

/* we took a full step. Thus add constraint nvl to the list of active */
/* constraints and update J and R */

	    ++(*nact);
	    iact[*nact-1] = nvl;

/* to update R we have to put the first nact-1 components of the d vector */
/* into column (nact) of R */

	    for (int i = 0; i < *nact - 1; i++) {
		work[iwrm + (*nact - 1) * *nact / 2 + i] = work[i];
/* L150: */
	    }

/* if now nact=n, then we just have to add the last element to the new */
/* row of R. */
/* Otherwise we use Givens transformations to turn the vector d(nact:n) */
/* into a multiple of the first unit vector. That multiple goes into the */
/* last element of the new row of R and J is accordingly updated by the */
/* Givens transformations. */

	    l = iwrm + (*nact) * (*nact + 1) / 2;
	    if (*nact == *n) {
		work[l-1] = work[*n-1];
	    } else {
		for (int i = *n - 1; i >= *nact; i--) {

/* we have to find the Givens rotation which will reduce the element */
/* (l1) of d to zero. */
/* if it is already zero we don't have to do anything, except of */
/* decreasing l1 */

		    if (work[i] == 0.) {
			goto L160;
		    }
/* Computing MAX */
		    d__3 = (d__1 = work[i - 1], abs(d__1)), d__4 = (d__2 = 
			    work[i], abs(d__2));
		    gc = max(d__3,d__4);
/* Computing MIN */
		    d__3 = (d__1 = work[i - 1], abs(d__1)), d__4 = (d__2 = 
			    work[i], abs(d__2));
		    gs = min(d__3,d__4);
		    d__1 = gc * sqrt((gs / gc) * (gs / gc) + 1);
		    temp = d_sign(&d__1, &work[i - 1]);
		    gc = work[i - 1] / temp;
		    gs = work[i] / temp;

/* The Givens rotation is done with the matrix (gc gs, gs -gc). */
/* If gc is one, then element (i) of d is zero compared with element */
/* (l1-1). Hence we don't have to do anything. */
/* If gc is zero, then we just have to switch column (i) and column (i-1) */
/* of J. Since we only switch columns in J, we have to be careful how we */
/* update d depending on the sign of gs. */
/* Otherwise we have to apply the Givens rotation to these columns. */
/* The i-1 element of d has to be updated to temp. */

		    if (gc == 1.) {
			goto L160;
		    }
		    if (gc == 0.) {
			work[i - 1] = gs * temp;
			for (int j = 0; j < *n; j++) {
			    temp = dmat[j + (i - 1) * dmat_dim1];
			    dmat[j + (i - 1) * dmat_dim1] = dmat[j + i * dmat_dim1];
			    dmat[j + i * dmat_dim1] = temp;
/* L170: */
			}
		    } else {
			work[i - 1] = temp;
			nu = gs / (gc + 1.);
			for (int j = 0; j < *n; j++) {
			    temp = gc * dmat[j + (i - 1) * dmat_dim1] + gs * dmat[j + i * dmat_dim1];
			    dmat[j + i * dmat_dim1] = nu * (dmat[j + (i - 1) * dmat_dim1] + temp) - dmat[j + i * dmat_dim1];
			    dmat[j + (i - 1) * dmat_dim1] = temp;
/* L180: */
			}
		    }
L160:
		    ;
		}

/* l is still pointing to element (nact,nact) of the matrix R. */
/* So store d(nact) in R(nact,nact) */
		work[l-1] = work[*nact-1];
	    }
	} else {

/* we took a partial step in dual space. Thus drop constraint it1, */
/* that is, we drop the it1-th active constraint. */
/* then we continue at step 2(a) (marked by label 55) */
/* but since the fit changed, we have to recalculate now "how much" */
/* the fit violates the chosen constraint now. */

	    sum = -bvec[nvl-1];
	    for (int j = 0; j < *n; j++) {
		sum += sol[j] * amat[j + (nvl-1) * amat_dim1];
/* L190: */
	    }
	    if (nvl > *meq) {
		work[iwsv + nvl-1] = sum;
	    } else {
		work[iwsv + nvl-1] = -abs(sum);
		if (sum > 0.) {
		    for (int j = 0; j < *n; j++) {
			amat[j + (nvl-1) * amat_dim1] = -amat[j + (nvl-1) * amat_dim1]
				;
/* L191: */
		    }
		    bvec[nvl-1] = -bvec[nvl-1];
		}
	    }
	    goto L700;
	}
    }
    goto L50;

/* Drop constraint it1 */

L700:

/* if it1 = nact it is only necessary to update the vector u and nact */

    if (it1 == *nact) {
	goto L799;
    }

/* After updating one row of R (column of J) we will also come back here */

L797:

/* we have to find the Givens rotation which will reduce the element */
/* (it1+1,it1+1) of R to zero. */
/* if it is already zero we don't have to do anything except of updating */
/* u, iact, and shifting column (it1+1) of R to column (it1) */
/* l  will point to element (1,it1+1) of R */
/* l1 will point to element (it1+1,it1+1) of R */

    l = iwrm + it1 * (it1 + 1) / 2 + 1;
    l1 = l + it1;
    if (work[l1-1] == 0.) {
	goto L798;
    }
/* Computing MAX */
    d__3 = (d__1 = work[l1-1 - 1], abs(d__1)), d__4 = (d__2 = work[l1-1], abs(
	    d__2));
    gc = max(d__3,d__4);
/* Computing MIN */
    d__3 = (d__1 = work[l1-1 - 1], abs(d__1)), d__4 = (d__2 = work[l1-1], abs(
	    d__2));
    gs = min(d__3,d__4);
    d__1 = gc * sqrt((gs / gc) * (gs / gc) + 1);
    temp = d_sign(&d__1, &work[l1-1 - 1]);
    gc = work[l1-1 - 1] / temp;
    gs = work[l1-1] / temp;

/* The Givens rotatin is done with the matrix (gc gs, gs -gc). */
/* If gc is one, then element (it1+1,it1+1) of R is zero compared with */
/* element (it1,it1+1). Hence we don't have to do anything. */
/* if gc is zero, then we just have to switch row (it1) and row (it1+1) */
/* of R and column (it1) and column (it1+1) of J. Since we swithc rows in */
/* R and columns in J, we can ignore the sign of gs. */
/* Otherwise we have to apply the Givens rotation to these rows/columns. */

    if (gc == 1.) {
	goto L798;
    }
    if (gc == 0.) {
	for (int i = it1 + 1; i <= *nact; i++) {
	    temp = work[l1-1 - 1];
	    work[l1-1 - 1] = work[l1-1];
	    work[l1-1] = temp;
	    l1 += i;
/* L710: */
	}
	for (int i = 0; i < *n; i++) {
	    temp = dmat[i + (it1-1) * dmat_dim1];
	    dmat[i + (it1-1) * dmat_dim1] = dmat[i + (it1) * dmat_dim1];
	    dmat[i + (it1) * dmat_dim1] = temp;
/* L711: */
	}
    } else {
	nu = gs / (gc + 1.);
	for (int i = it1 + 1; i <= *nact; i++) {
	    temp = gc * work[l1-1 - 1] + gs * work[l1-1];
	    work[l1-1] = nu * (work[l1-1 - 1] + temp) - work[l1-1];
	    work[l1-1 - 1] = temp;
	    l1 += i;
/* L720: */
	}
	for (int i = 0; i < *n; i++) {
	    temp = gc * dmat[i + (it1-1) * dmat_dim1] + gs * dmat[i + (it1) * dmat_dim1];
	    dmat[i + (it1) * dmat_dim1] = nu * (dmat[i + (it1-1) * dmat_dim1] + temp) - dmat[i + (it1) * dmat_dim1];
	    dmat[i + (it1-1) * dmat_dim1] = temp;
/* L721: */
	}
    }

/* shift column (it1+1) of R to column (it1) (that is, the first it1 */
/* elements). The posit1on of element (1,it1+1) of R was calculated above */
/* and stored in l. */

L798:
    for (int i = 0; i < it1; i++) {
	work[l-1 - it1 + i] = work[l-1 + i];
/* L730: */
    }

/* update vector u and iact as necessary */
/* Continue with updating the matrices J and R */

    work[iwuv + it1-1] = work[iwuv + it1];
    iact[it1-1] = iact[it1];
    ++it1;
    if (it1 < *nact) {
	goto L797;
    }
L799:
    work[iwuv + *nact-1] = work[iwuv + *nact];
    work[iwuv + *nact] = 0.;
    iact[*nact-1] = 0;
    --(*nact);
    ++(*pIterDeleted);
    goto L55;
L999:
    return 0;
} /* qpgen2_ */

