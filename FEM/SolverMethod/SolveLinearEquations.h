//
// Created by admin on 2023/1/14.
//

#ifndef MAIN_CPP_SOLVELINEAREQUATIONS_H
#define MAIN_CPP_SOLVELINEAREQUATIONS_H


#include <iostream>
#include <vector>
#include <valarray>

using namespace std;

class SolveLinearEquations {
protected:
    size_t n{};
    vector<vector<double>> A;
    valarray<double> b;

public:
    SolveLinearEquations() = default;

    SolveLinearEquations(vector<vector<double>> &A, valarray<double> &b, size_t n);

    ~SolveLinearEquations() = default;

};

#endif //MAIN_CPP_SOLVELINEAREQUATIONS_H
