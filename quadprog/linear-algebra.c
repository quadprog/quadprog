double sqrt(double);

void axpy(int n, double a, double x[], double y[])
{
    for (int i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

double dot(int n, double x[], double y[])
{
    double result = 0.0;
    for (int i =0; i < n; i++) {
        result += x[i] * y[i];
    }
    return result;
}

void scal(int n, double a, double x[])
{
    for (int i = 0; i < n; i++) {
        x[i] *= a;
    }
}

/*
 * Compute a * b, where a is upper triangular.
 * The result is written into b.
 */
void triangular_multiply(int n, double a[], double b[]) {
    for (int j = 0; j < n; j++) {
        axpy(j, b[j], &a[j * n], b);
        b[j] *= a[j + j * n];
    }
}

/*
 * Compute transpose(a) * b, where a is upper triangular.
 * The result is written into b.
 */
void triangular_multiply_transpose(int n, double a[], double b[]) {
    for (int j = n - 1; j >= 0; j--) {
        b[j] *= a[j + j * n];
        b[j] += dot(j, b, &a[j * n]);
    }
}

/*
 * Solve a * x = b, where a is upper triangular.
 * The solution is written into b.
 */
void triangular_solve(int n, double a[], double b[]) {
    for (int k = n - 1; k >= 0; k--) {
        b[k] /= a[k + k * n];
        axpy(k, -b[k], &a[k * n], b);
    }
}

/*
 * Solve transpose(a) * x = b, where a is upper triangular.
 * The solution is written into b.
 */
void triangular_solve_transpose(int n, double a[], double b[]) {
    for (int k = 0; k < n; k++) {
        b[k] -= dot(k, &a[k * n], b);
        b[k] /= a[k + k * n];
    }
}

/*
 * Invert a, where a is upper triangular.
 * The inverse is written into a.
 */
void triangular_invert(int n, double a[]) {
    for (int k = 0; k < n; k++) {
        a[k + k * n] = 1.0 / a[k + k * n];
        scal(k, -a[k + k * n], &a[k * n]);

        for (int j = k + 1; j < n; j++) {
            axpy(k, a[k + j * n], &a[k * n], &a[j * n]);
            a[k + j * n] *= a[k + k * n];
        }
    }
}

/*
 * Find the upper triangular matrix r such that a = transpose(r) * r, where a is positive definite.
 * The result is written into the upper triangle of a.
 * Returns: 0 if successful;
 *          j>0 if the leading jth minor of a is not positive definite.
 */
int cholesky(int n, double a[]) {
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < j; k++) {
            a[k + j * n] =
                (a[k + j * n] - dot(k, &a[k * n], &a[j * n])) / a[k + k * n];
        }

        double s = a[j + j * n] - dot(j, &a[j * n], &a[j * n]);

        if (s <= 0.0) {
            return j + 1;
        }

        a[j + j * n] = sqrt(s);
    }

    return 0;
}