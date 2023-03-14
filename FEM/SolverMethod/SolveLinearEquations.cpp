//
// Created by admin on 2023/1/14.
//

#include "SolveLinearEquations.h"

SolveLinearEquations::SolveLinearEquations(vector<vector<double>> &A, valarray<double> &b, size_t n) {
    this->A = A;
    this->b = b;
    this->n = n;
}


