#ifndef GAUSSEPP_CPP_MATRIX_H
#define GAUSSEPP_CPP_MATRIX_H
#pragma once

#include <iostream>
#include <vector>

using namespace std;

class Matrix {
public:
    size_t m_row;
    size_t m_col;
    vector<vector<double>> m_element;
public:

    Matrix(int m, int n, const vector<vector<double>> &A);

    ~Matrix() = default;

    Matrix(const Matrix &rhs) = default;

    Matrix &operator=(const Matrix &rhs) = default;

    double &element(int i, int j)  { return m_element[i][j]; }
    double element(int i, int j) const { return m_element[i][j]; }

};

#endif //GAUSSEPP_CPP_MATRIX_H
