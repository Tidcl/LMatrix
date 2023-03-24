//
// Created by admin on 2023/1/9.
//

#ifndef MAIN_CPP_DIRECTMETHOD_H
#define MAIN_CPP_DIRECTMETHOD_H

#include "../NMatrixDefine.h"
#include "SolveLinearEquations.h"

using namespace std;

class DirectMethod: public SolveLinearEquations{
public:
    DirectMethod(matrixv &A, valarray<double> &b, size_t n);

    ~DirectMethod() = default;

public:
    valarray<double> gauss_solver();

    valarray<double> tdma_solver();

};


#endif //MAIN_CPP_DIRECTMETHOD_H
