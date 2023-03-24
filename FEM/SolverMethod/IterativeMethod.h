//
// Created by admin on 2023/1/14.
//

#ifndef MAIN_CPP_ITERATIVEMETHOD_H
#define MAIN_CPP_ITERATIVEMETHOD_H

#include <iostream>
#include "../NMatrixDefine.h"
#include "SolveLinearEquations.h"

using namespace std;

class IterativeMethod: public SolveLinearEquations{
private:
    double eps;
    int max_iter;
public:
    IterativeMethod(matrixv &A, valarray<double> &b, size_t n,double eps,int max_iter);
public:
    valarray<double> gmres_solver();

    valarray<double> bicgstab_solver();
};


#endif //MAIN_CPP_ITERATIVEMETHOD_H
