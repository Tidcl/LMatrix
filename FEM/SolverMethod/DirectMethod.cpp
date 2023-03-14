//
// Created by admin on 2023/1/9.
//

#include "DirectMethod.h"

DirectMethod::DirectMethod(vector<vector<double>> &A, valarray<double> &b, size_t n) {
    this->A = A;
    this->b = b;
    this->n = n;
}


valarray<double> DirectMethod::gauss_solver() {
    //    pivoting
    valarray<double> x(n);
    int row;
    for (int k = 0; k < n - 1; ++k) {
        double maxnum = 0;
        for (int i = k; i < n; ++i) {
            if (fabs(A[i][k]) > maxnum) {
                maxnum = fabs(A[i][k]);
                row = i;
            }
        }
        if (k != row) {
            for (int i = 0; i < n; ++i) {
//                mySwap(A[row][i], A[k][i]);
                swap(A[row][i], A[k][i]);
            }
//            mySwap(b[row], b[k]);
            swap(b[row], b[k]);
        }
    }
//    eliminate
    double m;
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            m = A[i][k] / A[k][k];
            for (int j = 0; j < n; ++j) {
                A[i][j] -= m * A[k][j];
            }
            b[i] -= m * b[k];
        }
    }
//    backsubstitute
    for (int k = int(n) - 1; k >= 0; --k) {
        double sum = 0;
        for (int j = k; j <= n - 1; ++j) {
            sum += A[k][j] * x[j];
        }
        x[k] = (b[k] - sum) / A[k][k];
    }
    return x;
}

valarray<double> DirectMethod::tdma_solver() {
    valarray<double> coe_a(n), coe_b(n), coe_c(n);
    valarray<double> coe_d = b;
    for (int i = 0; i < n; ++i) coe_b[i] = A[i][i];
    for (int i = 1; i < n; ++i) coe_a[i] = A[i][i - 1];
    for (int i = 0; i < n - 1; ++i) coe_c[i] = A[i][i + 1];

    valarray<double> x(n);

    for (int i = 0; i < n; ++i) {
        if (coe_b[i] + 1.0 == 1.0) {
            cout << "Matrix singular, exit" << endl;
            exit(1);
        }
    }

    coe_c[0] = coe_c[0] / coe_b[0];
    coe_d[0] = coe_d[0] / coe_b[0];

    for (int i = 1; i < n; ++i) {
        double m = 1.0 / (coe_b[i] - coe_a[i] * coe_c[i - 1]);
        coe_c[i] = m * coe_c[i];
        coe_d[i] = m * (coe_d[i] - coe_a[i] * coe_d[i - 1]);
    }

    x[n - 1] = coe_d[n - 1];
    for (int i = int(n) - 2; i >= 0; --i) x[i] = coe_d[i] - coe_c[i] * x[i + 1];

    return x;
}