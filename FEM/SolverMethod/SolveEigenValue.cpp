//
// Created by admin on 2023/1/10.
//

#include "SolveEigenValue.h"
#include "../Tools/PrintTools.h"


SolveEigenValue::SolveEigenValue(vector<vector<double>> &A, int &max_iter, double &eps) {
    this->A = A;
    this->max_iter = max_iter;
    this->eps = eps;
}

void SolveEigenValue::power_method() {
    int n = int(A.size());
    vector<double> m_x(n, 1.0);
    vector<double> m_y(n, 1.0);

    vector<double> y_next(n, 1.0);
    int iterateNum = 0;
    double sum;
    double maxEigen;

    do {
        iterateNum++;
        for (int i = 0; i < n; ++i) m_y[i] = y_next[i];
        for (int i = 0; i < n; ++i) {
            sum = 0;
            for (int j = 0; j < n; ++j) sum += this->A[i][j] * m_y[j];
            m_x[i] = sum;
        }

        maxEigen = 0.0;
        for (int i = 0; i < n; ++i) {
            if (fabs(maxEigen) < fabs(m_x[i])) maxEigen = m_x[i];
        }
        for (int i = 0; i < n; ++i) y_next[i] = m_x[i] / maxEigen;

        sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += fabs(m_y[i] - y_next[i]);
        }

        if (iterateNum < max_iter) {
            if (fabs(sum) < eps) {
                cout << "Tolerance = " << sum << endl;
                break;
            }
        } else {
            cout << "Iterate Num =" << iterateNum << endl;
            break;
        }
    } while (true);

    int accuracy = 8;
    MyPrint myPrint(accuracy);

    cout << "Max Eigen Value = " << maxEigen << endl;
    cout << "Eigen Vector:" << '\n';
    myPrint.print_vector(m_y);
}

void SolveEigenValue::jacobi_rush() {
    int i, j, p, q, u, w, t, s;
    double ff, fm, cn, sn, omega, x, y, d;

    valarray<double> Aa = vtr2d_to_valary1d(A);
    int n = int(sqrt(double(Aa.size())));

    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            if (Aa[i * n + j] - Aa[j * n + i] >= 1e-5) {
                cout << "Matrix is not match!" << endl;
            }
        }
    }

    valarray<double> v(Aa.size());
    for (i = 0; i < n; ++i) {
        v[i * n + i] = 1.0;
        for (j = 0; j < n; ++j) {
            if (i != j) {
                v[i * n + j] = 0.0;
            }
        }
    }

    ff = 0.0;
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            d = Aa[i * n + j];
            ff += pow(d, 2);
        }
    }
    ff = sqrt(2.0 * ff);
    ff /= (1.0 * double(n));
    while (ff >= this->eps) {
        d = 0.0;
        for (i = 1; i < n and d <= ff; ++i) {
            for (j = 0; j < i and d <= ff; ++j) {
                d = fabs(Aa[i * n + j]);
                p = i, q = j;
            }
        }
        if (d <= ff) {
            ff /= (1.0 * double(n));
        } else {
            u = p * int(n) + q;
            w = p * int(n) + p;
            t = q * int(n) + p;
            s = q * int(n) + q;
            x = -Aa[u];
            y = (Aa[s] - Aa[w]) / 2.0;
            omega = x / sqrt(pow(x, 2) + pow(y, 2));
            if (y < 0.0) {
                omega = -omega;
            }
            sn = 1.0 + sqrt(1.0 - pow(omega, 2));
            sn = omega / sqrt(2.0 * sn);
            cn = sqrt(1.0 - sn * sn);
            fm = Aa[w];
            Aa[w] = fm * pow(cn, 2) + Aa[s] * pow(sn, 2) + Aa[u] * omega;
            Aa[s] = fm * pow(sn, 2) + Aa[s] * pow(cn, 2) - Aa[u] * omega;
            Aa[u] = 0.0;
            Aa[t] = 0.0;
            for (j = 0; j < n; ++j) {
                if (j != p and j != q) {
                    u = p * int(n) + j;
                    w = q * int(n) + j;
                    fm = Aa[u];
                    Aa[u] = fm * cn + Aa[w] * sn;
                    Aa[w] = -fm * sn + Aa[w] * cn;
                }
            }
            for (i = 0; i < n; ++i) {
                if (i != p and i != q) {
                    u = i * int(n) + p;
                    w = i * int(n) + q;
                    fm = Aa[u];
                    Aa[u] = fm * cn + Aa[w] * sn;
                    Aa[w] = -fm * sn + Aa[w] * cn;
                }
            }
            for (i = 0; i < n; ++i) {
                u = i * int(n) + p;
                w = i * int(n) + q;
                fm = v[u];
                v[u] = fm * cn + v[w] * sn;
                v[w] = -fm * sn + v[w] * cn;
            }
        }
    }

    vector<vector<double>> AA = valary1d_to_vtr2d(Aa);
    vector<vector<double>> VV = valary1d_to_vtr2d(v);

    int accuracy = 8;
    MyPrint myPrint(accuracy);
    cout << "Eigen Value Matrix" << endl;
    myPrint.print_matrix(AA);
    cout << "Eigen Vector Matrix" << endl;
    myPrint.print_matrix(VV);
}