//
// Created by admin on 2023/1/10.
//

#ifndef MAIN_CPP_SOLVEEIGENVALUE_H
#define MAIN_CPP_SOLVEEIGENVALUE_H

#include <iostream>
#include "../NMatrixDefine.h"
#include "../Tools/TypeConversion.h"

using namespace std;


class SolveEigenValue {
private:
    vector<vector<double>> A;
    int max_iter;
    double eps;
public:
    SolveEigenValue(vector<vector<double>> &A, int &max_iter, double &eps);

    ~SolveEigenValue() = default;

public:
    void power_method();

    void jacobi_rush();
};


#endif //MAIN_CPP_SOLVEEIGENVALUE_H
