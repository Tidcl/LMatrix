#include "Gmres.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template<class Tv>
Gmres<Tv>::Gmres() = default;

template<class Tv>
Gmres<Tv>::~Gmres() = default;

template<class Tv>
Gmres<Tv>::Gmres(const int n) {
    m_size = n;
    m_b.resize(n);
    m_mtx.resize(n);
    for (size_t i = 0; i < n; ++i) {
        m_mtx[i].resize(n);
    }
}

template<class Tv>
Gmres<Tv>::Gmres(const int n, const int m, const vector<Tv> &A, const Tv &b) {
    m_size = n;
    m_b = b;
    m_mtx = A;
}

template<class Tv>
void Gmres<Tv>::printMatrix(vector<Tv> &xy) {
    for (int i = 0; i < xy.size(); ++i) {
        for (int j = 0; j < xy[0].size(); ++j) {
            cout << xy[i][j] << "\t";
        }
        cout << endl;
    }
}

template<class Tv>
void Gmres<Tv>::printVector(Tv &x) {
    for (int i = 0; i < x.size(); ++i) {
        cout << x[i] << "\t";
    }
    cout << endl;
}

template<class Tv>
double dot(const Tv &x, const Tv &y) {
    double sum = 0;
    if (x.size() != y.size()) {
        cout << "size not match\n";
        exit(1);
    } else {
        for (int i = 0; i < x.size(); ++i) {
            sum += x[i] * y[i];
        }
    }
    return sum;
}

template<class Tv>
double norm2(const Tv &a) {
    return sqrt(dot(a, a));
}

template<class Tv>
Tv operator*(const double C, const Tv &x) {
    Tv res(x.size());
    for (int i = 0; i < x.size(); ++i) res[i] = C * x[i];
    return res;

}

template<class Tv>
Tv operator*(const vector<Tv> &A, const Tv &x) {
    Tv tm(x.size());
    if (A[0].size() != x.size()) {
        cout << "size not match\n";
        exit(1);
    } else {
        for (int i = 0; i < A.size(); ++i) {
            for (int j = 0; j < A[0].size(); ++j) {
                tm[i] = tm[i] + A[i][j] * x[j];
            }
        }
    }
    return tm;
}

template<class Tv>
Tv operator-(const Tv &x, const Tv &y) {
    Tv res(x.size());
    if (x.size() != y.size()) {
        cout << "size not match\n";
        exit(1);
    } else {
        for (int i = 0; i < x.size(); ++i) res[i] = x[i] - y[i];
    }
    return res;
}


template<class Tv>
void Gmres<Tv>::Arnoldi(const int m, vector<Tv> &h, vector<Tv> &v, double beta) {
    size_t n = this->m_size;
    Tv x(n);
    Tv r(n);
//    printMatrix(m_mtx);
//    printVector(x);
    Tv Ar = this->m_mtx * x;
//    printVector(Ar);
    r = this->m_b - Ar;
//    printVector(r);
    Tv w(n);
    v.push_back((1.0 / beta) * r);
    for (int k = 0; k < m; ++k) {
        w = this->m_mtx * v[k];
        for (int i = 0; i <= k; ++i) {
            h[k][i] = dot(w, v[i]);
            w = w - h[k][i] * v[i];
        }
        h[k][k + 1] = norm2(w);
        if (h[k][k + 1] + 1.0 == 1.0) {
            cout << "h[k][k+1] = \n";
            exit(1);
        }
        v.push_back((1.0 / h[k][k + 1]) * w);
    }
//    printMatrix(v);
}

template<class Tv>
void Gmres<Tv>::Givens(const int m, vector<Tv> &h, Tv &g, double beta, Tv &x, vector<Tv> &v) {
    double eps = 1.0e-10;
    double stp = norm2(m_b) * eps;
    int iter = 1;
    int itmax = 6;

    while (iter <= itmax) {
        Tv cs(m), sn(m);
        g[0] = beta;
        int k;
        for (k = 0; k < m and iter <= itmax; ++k, ++iter) {
            for (int i = 0; i < k; ++i) {
                double temp = cs[i] * h[k][i] + sn[i] * h[k][i + 1];
                h[k][i + 1] = -sn[i] * h[k][i] + cs[i] * h[k][i + 1];
                h[k][i] = temp;
            }
            if (h[k][k + 1] + 1.0 == 1.0) {
                cs[k] = 1;
                sn[k] = 0;
            } else {
                double dothh = h[k][k] * h[k][k] + h[k][k + 1] * h[k][k + 1];
                double tpm = sqrt(dothh);
                cs[k] = h[k][k] / tpm;
                sn[k] = h[k][k + 1] / tpm;
            }
            double dothh1 = cs[k] * h[k][k] + sn[k] * h[k][k + 1];
            h[k][k] = dothh1;
            double tmpgk = cs[k] * g[k] + sn[k] * g[k + 1];
            g[k + 1] = -sn[k] * g[k] + cs[k] * g[k + 1];
            g[k] = tmpgk;

            if (fabs(g[k + 1]) <= stp) {
                k++;
                break;
            }
        }
        for (int i = k - 1; i >= 0; --i) {
            for (int j = i + 1; j < k; ++j) {
                g[i] -= h[j][i] * g[j];
            }
            g[i] /= h[i][i];
        }

        for (int i = 0; i < k; ++i) {
            Tv tm = g[i] * v[i];
            for (int j = 0; j < v[0].size(); ++j) x[j] = x[j] + tm[j];
        }
        beta = norm2(this->m_b - this->m_mtx * x);
        if (fabs(beta) <= stp) break;
    }
    printVector(x);
}

template<class Tv>
bool Gmres<Tv>::solve(const int m) {
    double eps = 1.0e-10;
    double stp = norm2(m_b) * eps;
    size_t n = this->m_size;
    Tv x(n);
    vector<vector<double>> v;
    vector<vector<double>> h(m, vector<double>(m + 1));
    Tv r(n);
    Tv g(m + 1);
    Tv Ar = this->m_mtx * x;
    r = this->m_b - Ar;
    double beta = norm2(r);

    Arnoldi(m, h, v, beta);
    Givens(m, h, g, beta, x, v);

    double epsilon = norm2(m_b - m_mtx * x);
    cout << "Residula in Gmres=" << epsilon << endl;
    if (eps < stp) {
        return true;
    } else {
        return false;
    }
}


int main() {
    cout.setf(ios::fixed);
    vector<double> b0 = {11, 9, 9, 9, 13, 17};
    vector<vector<double>> A0 = {{6, 5, 4, 3, 2, 1},
                                 {5, 6, 5, 4, 3, 2},
                                 {4, 5, 6, 5, 4, 3},
                                 {3, 4, 5, 6, 5, 4},
                                 {2, 3, 4, 5, 6, 5},
                                 {1, 2, 3, 4, 5, 6}};

    const int m = b0.size();
    Gmres<vector<double>> g1(6, 6, A0, b0);
//    g1.printVector(b0);
//    g1.printMatrix(A0);
    g1.solve(m);
    return 0;
}
