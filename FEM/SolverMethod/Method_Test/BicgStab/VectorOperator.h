#ifndef GAUSSEPP_CPP_VECTOROPERATOR_H
#define GAUSSEPP_CPP_VECTOROPERATOR_H
#pragma once

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


//add by tangxc 不建议在不满足的情况下打印错误或者直接退出，
//add by tangxc 作为工具类最好是返回一个空结果代表运算失败，打印等业务信息让用户自己去做
template<class Tv = vector<double>>
Tv operator*(const double C, const Tv &x) {
    Tv res(x.size());
    for (int i = 0; i < x.size(); ++i) res[i] = C * x[i];
    return res;
}

template<class Tv = vector<double>>
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

template<class Tv = vector<double>>
Tv operator-(const Tv &y, const Tv &x) {
    Tv res(x.size());
    if (x.size() != y.size()) {
        cout << "Size not match \n";
        exit(1);
    } else {
        for (int i = 0; i < x.size(); ++i) res[i] = y[i] - x[i];
        return res;
    }
}

template<class Tv = vector<double>>
double dot(const Tv &x, const Tv &y) {
    double sum = 0;
    if (x.size() != y.size()) {
        cout << "Size not match \n";
        exit(1);
    } else {
        for (int i = 0; i < x.size(); ++i) {
            sum += x[i] * y[i];
        }
    }
    return sum;
}

template<class Tv = vector<double>>
double norm2(const Tv &a){
    return sqrt(dot(a,a));
}

#endif //GAUSSEPP_CPP_VECTOROPERATOR_H
