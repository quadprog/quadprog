/* Table of constant values */

static int c__1 = 1;

/* Subroutine */ int dpori_(double *a, int *lda, int *n)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int j, k;
    double t;
    int kp1;
    extern /* Subroutine */ int dscal_(int *, double *, double *,
	    int *), daxpy_(int *, double *, double *, int
	    *, double *, int *);


/*     dpori computes the inverse of the factor of a */
/*     double precision symmetric positive definite matrix */
/*     using the factors computed by dpofa. */

/*     modification of dpodi by BaT 05/11/95 */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the output  a  from dpofa */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       if dpofa was used to factor  a  then */
/*                dpodi produces the upper half of inverse(a) . */
/*                elements of  a  below the diagonal are unchanged. */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if dpoco or dpofa has set info .eq. 0 . */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */
/*     modified by Berwin A. Turlach 05/11/95 */

/*     subroutines and functions */

/*     blas daxpy,dscal */
/*     fortran mod */

/*     internal variables */


/*     compute inverse(r) */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
	t = -a[k + k * a_dim1];
	i__2 = k - 1;
	dscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1];
	    a[k + j * a_dim1] = 0.;
	    daxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }
    return 0;
} /* dpori_ */

/* Subroutine */ int dposl_(double *a, int *lda, int *n,
	double *b)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int k;
    double t;
    int kb;
    extern double ddot_(int *, double *, int *, double *,
	    int *);
    extern /* Subroutine */ int daxpy_(int *, double *, double *,
	    int *, double *, int *);


/*     dposl solves the double precision symmetric positive definite */
/*     system a * x = b */
/*     using the factors computed by dpoco or dpofa. */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the output from dpoco or dpofa. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        b       double precision(n) */
/*                the right hand side vector. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal.  technically this indicates */
/*        singularity but it is usually caused by improper subroutine */
/*        arguments.  it will not occur if the subroutines are called */
/*        correctly and  info .eq. 0 . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call dpoco(a,lda,n,rcond,z,info) */
/*           if (rcond is too small .or. info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call dposl(a,lda,n,c(1,j)) */
/*        10 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */

/*     internal variables */


/*     solve trans(r)*y = b */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = ddot_(&i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	b[k] = (b[k] - t) / a[k + k * a_dim1];
/* L10: */
    }

/*     solve r*x = y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* dposl_ */

