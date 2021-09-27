double hypot(double, double);

/*
 * Apply orthogonal transformations to a to bring the components beyond the rth to zero.
 * Append the result to R as a final column.
 * Apply the same orthogonal transformations to the columns of Q.
 *
 * Note that R is an (r-1) by (r-1) upper triangular matrix stored as packed columns.
 * So the input size is (r-1)*r/2 and the output size is r*(r+1)/2.
 */
void qr_insert(int n, int r, double a[], double Q[], double R[]) {
    for (int i = n - 1; i >= r; i--) {
        // On this iteration, reduce a[i] to zero.

        if (a[i] == 0.) {
            continue;
        }

        if (a[i - 1] == 0.) {
            // Simply swap

            a[i - 1] = a[i];
            for (int j = 0; j < n; j++) {
                double temp = Q[j + (i - 1) * n];
                Q[j + (i - 1) * n] = Q[j + i * n];
                Q[j + i * n] = temp;
            }
        } else {
            // Compute a Givens rotation.

            double h = hypot(a[i - 1], a[i]);

            if (a[i - 1] < 0.) {
                h = -h;
            }

            double gc = a[i - 1] / h;
            double gs = a[i] / h;
            double nu = a[i] / (a[i - 1] + h); // this saves a fourth multiplication in the inner loop below

            a[i - 1] = h;

            for (int j = 0; j < n; j++) {
                double temp = gc * Q[j + (i - 1) * n] + gs * Q[j + i * n];
                Q[j + i * n] = nu * (Q[j + (i - 1) * n] + temp) - Q[j + i * n];
                Q[j + (i - 1) * n] = temp;
            }
        }
    }

    for (int i = 0; i < r; i++) {
        R[(r - 1) * r / 2 + i] = a[i];
    }
}

/*
 * Drop the col-th column of R.
 * Apply orthogonal transformations to the rows of R to restore R to upper triangular form.
 * Apply the same orthogonal transformations to the columns of Q.
 *
 * Note that R is an r by r upper triangular matrix stored as packed columns.
 * So the input size is r*(r+1)/2 and the output size is (r-1)*r/2.
 */
void qr_delete(int n, int r, int col, double Q[], double R[]) {
    for (int i = col; i < r; i++) {
        // On this iteration, reduce the (i+1, i+1) element of R to zero,
        // and then move column (i+1) to position i.

        // R[l] is the (i+1, i+1) element of R.
        int l = (i + 1) * (i + 2) / 2 - 1;

        if (R[l] == 0.) {
            continue;
        }

        if (R[l - 1] == 0.) {
            // Simply swap.

            for (int j = i + 1; j <= r; j++) {
                double temp = R[l - 1];
                R[l - 1] = R[l];
                R[l] = temp;
                l += j;
            }

            for (int j = 0; j < n; j++) {
                double temp = Q[j + (i - 1) * n];
                Q[j + (i - 1) * n] = Q[j + i * n];
                Q[j + i * n] = temp;
            }
        } else {
            // Compute a Givens rotation.

            double h = hypot(R[l - 1], R[l]);

            if (R[l - 1] < 0.0) {
                h = -h;
            }

            double gc = R[l - 1] / h;
            double gs = R[l] / h;
            double nu = R[l] / (R[l - 1] + h); // this saves a fourth multiplication in the inner loop below

            for (int j = i + 1; j <= r; j++) {
                double temp = gc * R[l - 1] + gs * R[l];
                R[l] = nu * (R[l - 1] + temp) - R[l];
                R[l - 1] = temp;
                l += j;
            }

            for (int j = 0; j < n; j++) {
                double temp = gc * Q[j + (i - 1) * n] + gs * Q[j + i * n];
                Q[j + i * n] = nu * (Q[j + (i - 1) * n] + temp) - Q[j + i * n];
                Q[j + (i - 1) * n] = temp;
            }
        }

        for (int j = 0; j < i; j++) {
            R[(i - 1) * i / 2 + j] = R[i * (i + 1) / 2 + j];
        }
    }
}