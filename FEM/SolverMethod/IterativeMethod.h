//
// Created by admin on 2023/1/14.
//

#ifndef MAIN_CPP_ITERATIVEMETHOD_H
#define MAIN_CPP_ITERATIVEMETHOD_H

#include <iostream>
#include <vector>
#include <valarray>
#include "SolveLinearEquations.h"

using namespace std;

class IterativeMethod: public SolveLinearEquations{
private:
    size_t n{};
    vector<vector<double>> A;
    valarray<double> b;
    double eps;
    int max_iter;
public:
    IterativeMethod(vector<vector<double>> &A, valarray<double> &b, size_t n,double eps,int max_iter);
public:
    valarray<double> gmres_solver();

    valarray<double> bicgstab_solver();
};


#endif //MAIN_CPP_ITERATIVEMETHOD_H
