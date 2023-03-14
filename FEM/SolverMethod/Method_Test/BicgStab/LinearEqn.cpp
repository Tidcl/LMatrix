#include "LinearEqn.h"
#include "VectorOperator.h"

using namespace std;

template<class Tv = vector<double>>
Tv operator*(const Matrix &A, const Tv &x) {
    Tv tm(x.size());
    if (A.m_row != x.size()) {
        cout << "Size not match \n";
        exit(1);
    } else {
        for (int i = 0; i < A.m_row; ++i) {
            for (int j = 0; j < A.m_col; ++j) {
                tm[i] = tm[i] + A.element(i, j) * x[j];
            }
        }
    }
    return tm;
}

void LinearEqn::printVector(vector<double> &v) {
    cout << "Vector is: \n";
    for (size_t i = 0; i < v.size(); ++i) {
        cout << v[i] << endl;
    }
    cout << endl;
}

void LinearEqn::printMatrix(vector<vector<double>> &A) {
    cout << "----------Matrix-----------" << endl;
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            cout << m_element[i][j] << '\t';
        }
        cout << endl;
    }
}


void LinearEqn::printSolution() {
    cout << "Solution is: \n";
    for (size_t i = 0; i < m_solution.size(); ++i) {
        cout << m_solution[i] << endl;
    }
    cout << endl;
}

LinearEqn::LinearEqn(int m, int n, const vector<vector<double>> &A, const vector<double> &b) : Matrix(m, n, A) {
    m_b = b;
//    cout<<m<<endl;
//    cout<<n<<endl;
}

bool LinearEqn::BicgStab(const double tol, const int maxIter) {
//    printMatrix(m_element);
//    printVector(m_b);

    size_t bsize = m_b.size();
    vector<double> x(bsize, 0), p(bsize), v(bsize), s(bsize), t(bsize);
    vector<double> r = m_b - (*this) * x;
    vector<double> r_bar(r);

    double rho = 1, alpha = 1, w = 1;
    double rho_new, beta;
    double residual;

    int num = 1;
    while (num < maxIter) {
        rho_new = dot(r_bar, r);
        if (rho_new + 1.0 == 1.0) break;
        if (num == 1) {
            p = r;
        } else {
            beta = (rho_new / rho) * (alpha / w);
            p = r + beta * (p - w * v);
        }

        v = (*this) * p;
        alpha = rho_new / (dot(r_bar, v));
        s = r - alpha * v;

        if (norm2(s) < tol) {
            x = x + alpha * p;
            break;
        }

        t = (*this) * s;
        w = dot(t, s) / dot(t, t);
        x = x + alpha * p + w * s;
        r = s - w * t;

        residual = norm2(s) / norm2(m_b);
        if ((residual < tol) or (w + 1.0 == 1.0)) { break; }

        num++;
        rho = rho_new;
    }
    m_solution = x;
    cout << "iter num = " << num << endl;
    cout << scientific << "residual = " << residual << endl;
    cout << fixed;
    if (residual > tol) {
        cout << "can not converge \n";
        return false;
    } else {
        return true;
    }
}
