#include "Matrix.h"

Matrix::Matrix(int m, int n, const vector<vector<double>> &A) {
    this->m_row = m;
    this->m_col = n;
    this->m_element = A;
}