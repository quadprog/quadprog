#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif

extern void F_FUNC(qpgen2,QPGEN2)(double*,double*,int*,int*,double*,double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,double*,int*);

void qpgen2(double* dmat, double* dvec,int* fddmat,int* n,
            double* sol, double* lagr, double* crval,
            double* amat, double* bvec, int* fdamat,
            int* q,int* meq, int* iact, int* nact, int* iter,
            double* work, int* ierr) {
    F_FUNC(qpgen2,QPGEN2)(dmat, dvec, fddmat, n, sol, lagr, crval,
        amat, bvec, fdamat, q, meq, iact, nact, iter, work, ierr);
}
