#ifndef GAUSSEPP_CPP_LINEAREQN_H
#define GAUSSEPP_CPP_LINEAREQN_H
#pragma once

#include <iostream>
#include "Matrix.h"
#include <vector>

using namespace std;


class LinearEqn : public Matrix {
public:
    LinearEqn(int m, int n, const vector<vector<double>> &A, const vector<double> &b);
    ~LinearEqn() = default;

public:
    vector<double> m_b;
    vector<double> m_solution;
public:
    void printSolution();
    static void printVector(vector<double> &b);
    void printMatrix(vector<vector<double>> &A);
    bool BicgStab(double tol, int maxIter);
};

#endif //GAUSSEPP_CPP_LINEAREQN_H
