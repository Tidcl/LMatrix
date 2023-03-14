//
// Created by admin on 2023/1/14.
//

#include "../Tools/PrintTools.h"
#include "../Tools/TypeConversion.h"
#include "../Tools/SymbolsOperat.h"
#include "IterativeMethod.h"

IterativeMethod::IterativeMethod(vector<vector<double>> &A, valarray<double> &b, size_t n, double eps, int max_iter) {
    this->A = A;
    this->b = b;
    this->n = n;
    this->eps = eps;
    this->max_iter = max_iter;
}

valarray<double> IterativeMethod::gmres_solver() {


    double stp = norm2(b) * eps;
    int m = int(b.size());

    vector<double> x(n);
    vector<vector<double>> v;
    vector<vector<double>> h(m, vector<double>(m + 1));
    vector<double> r(n);
    vector<double> g(m + 1);
    vector<double> Ar = this->A * x;
    r = this->b - Ar;
    double beta = norm2(r);
    vector<double> w(n);

    v.push_back((1.0 / beta) * r);


    for (int k = 0; k < m; ++k) {
        w = this->A * v[k];
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


    int iter = 1;
    while (iter < max_iter) {
        vector<double> cs(m);
        vector<double> sn(m);
        g[0] = beta;
        int k;
        for (k = 0; k < m and iter < max_iter; ++k, ++iter) {
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
            vector<double> tm = g[i] * v[i];
            for (int j = 0; j < v[0].size(); ++j) x[j] = x[j] + tm[j];
        }
        beta = norm2(this->b - this->A * x);
        if (fabs(beta) <= stp) break;
    }

    cout << "Gmres Iter Num = " << iter << endl;

    return vtr1d_to_valary1d(x);
}


valarray<double> IterativeMethod::bicgstab_solver() {
    size_t bsize = b.size();
    vector<double> x(bsize, 0), p(bsize), v(bsize), s(bsize), t(bsize);
    vector<double> r = b - (this->A * x);
    vector<double> r_bar(r);

    double rho = 1, alpha = 1, w = 1;
    double rho_new, beta;
    double residual;

    int num = 1;
    while (num < max_iter) {
        rho_new = dot(r_bar, r);
        if (rho_new + 1.0 == 1.0) break;
        if (num == 1) {
            p = r;
        } else {
            beta = (rho_new / rho) * (alpha / w);
            p = r + beta * (p - w * v);
        }

        v = this->A * p;
        alpha = rho_new / (dot(r_bar, v));
        s = r - alpha * v;

        if (norm2(s) < eps) {
            x = x + alpha * p;
            break;
        }

        t = this->A * s;
        w = dot(t, s) / dot(t, t);
        x = x + alpha * p + w * s;
        r = s - w * t;

        residual = norm2(s) / norm2(b);
        if ((residual < eps) or (w + 1.0 == 1.0)) { break; }

        num++;
        rho = rho_new;
    }

    cout << "BICG_Stab Iter Num = " << num << endl;
//    cout << scientific << "residual = " << residual << endl;
//    cout << fixed;
//    if (residual > eps) {
//        cout << "residual > eps \n";
//        cout << "can not converge \n";
//    } else {
//        cout << "residual < eps \n";
//    }
    return vtr1d_to_valary1d(x);
}



