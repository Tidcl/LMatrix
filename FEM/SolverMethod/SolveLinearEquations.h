//
// Created by admin on 2023/1/14.
//

#ifndef MAIN_CPP_SOLVELINEAREQUATIONS_H
#define MAIN_CPP_SOLVELINEAREQUATIONS_H


#include <iostream>
#include "../NMatrixDefine.h"

using namespace std;

class SolveLinearEquations {
protected:
    size_t n{};
    matrixv A;
    valarray<double> b;

public:
    SolveLinearEquations() = default;

    SolveLinearEquations(matrixv &A, valarray<double> &b, size_t n);

    virtual ~SolveLinearEquations() = default;

};

#endif //MAIN_CPP_SOLVELINEAREQUATIONS_H
