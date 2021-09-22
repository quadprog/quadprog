double ddot_(int *n, double *dx, int *incx, double *dy,
	int *incy)
{
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    int i__, m, ix, iy, mp1;
    double dtemp;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     DDOT forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

	m = *n % 5;
	if (m != 0) {
	    i__1 = m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dtemp += dx[i__] * dy[i__];
	    }
	    if (*n < 5) {
		ret_val = dtemp;
		return ret_val;
	    }
	}
	mp1 = m + 1;
	i__1 = *n;
	for (i__ = mp1; i__ <= i__1; i__ += 5) {
	    dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + 
		    dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + 
		    dx[i__ + 4] * dy[i__ + 4];
	}
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
	    ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
	    iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dtemp += dx[ix] * dy[iy];
	    ix += *incx;
	    iy += *incy;
	}
    }
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

