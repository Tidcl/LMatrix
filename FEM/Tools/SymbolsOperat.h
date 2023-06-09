//
// Created by Administrater on 2023/1/8.
//

#ifndef RUNSCRIPY_PY_SYMBOLSOPERAT_H
#define RUNSCRIPY_PY_SYMBOLSOPERAT_H

#include "../NMatrixDefine.h"

using namespace std;

//add by tangxc 不建议在不满足的情况下打印错误或者直接退出，
//add by tangxc 作为工具类最好是返回一个空结果代表运算失败，打印等业务信息让用户自己去做
template<class Tv = vecd>
Tv operator*(const double C, const Tv &x) {
    Tv res(x.size(), 0);
    for (int i = 0; i < x.size(); ++i) res[i] = C * x[i];
    return res;
}


template<class Tv = vecd>
Tv operator+(const Tv &x, const double C) {
    Tv result(x.size(), 0);
    for (int i = 0; i < x.size(); ++i) result[i] = C + x[i];
    return result;
}

template<class Tv = vecd>
Tv operator+(const Tv &y, const Tv &x) {
    Tv result(x.size(), 0);
    if (x.size() != y.size()) {
        cout << "Vector size not match \n";
        exit(1);
    } else {
        for (int i = 0; i < x.size(); ++i) result[i] = y[i] + x[i];
        return result;
    }
}

template<class Tv = vecd>
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

template<class Tv = vecd>
vector<Tv> operator+(matrixv &A,matrixv &B) {
    matrixv tm(A.size(),vecd(A[0].size(),0));
    if (A[0].size() != B[0].size()) {
        cout << "size not match\n";
        exit(1);
    } else {
        for (int i = 0; i < A.size(); ++i) {
            for (int j = 0; j < A[0].size(); ++j) {
                tm[i][j] = tm[i][j] + (A[i][j] + B[i][j]);
            }
        }
    }
    return tm;
}

template<class Tv = vecd>
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

template<class Tv = vecd>
Tv operator-(const valarray<double> &x, const Tv &y) {
    Tv res(x.size());
    if (x.size() != y.size()) {
        cout << "size not match\n";
        exit(1);
    } else {
        for (int i = 0; i < x.size(); ++i) res[i] = x[i] - y[i];
    }
    return res;
}

template<class Tv= vecd>
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

template<class Tv= vecd>
double norm2(const Tv &a) {
    return sqrt(dot(a, a));
}


//template<class Tv = vecd>
//Tv operator*(const double C, const Tv &x) {
//    Tv res(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) res[i] = C * x[i];
//    return res;
//}
//
//template<class Tv = vecd>
//Tv operator+(const Tv &x, const double C) {
//    Tv result(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) result[i] = C + x[i];
//    return result;
//}

#endif //RUNSCRIPY_PY_SYMBOLSOPERAT_H
