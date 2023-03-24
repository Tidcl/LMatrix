//
// Created by admin on 2023/1/9.
//

#ifndef MAIN_CPP_TYPECONVERSION_H
#define MAIN_CPP_TYPECONVERSION_H

#include "../NMatrixDefine.h"

using namespace std;

static valarray<double> vtr2d_to_valary1d(matrixv &x2d) {
    size_t row_num = x2d.size();
    size_t col_num = x2d[0].size();
    valarray<double> x1d(row_num * col_num);
    for (int i = 0; i < row_num; ++i) {
        for (int j = 0; j < col_num; ++j) {
            x1d[row_num * j + i] = x2d[i][j];
        }
    }
    return x1d;
}


static matrixv valary1d_to_vtr2d(valarray<double> &x1d) {
    size_t tol_num = x1d.size();
    auto row_num = size_t(sqrt(double_t (tol_num)));
    size_t col_num = row_num;
    matrixv x2d(row_num, vecd(col_num, 0));
    for (int i = 0; i < row_num; ++i) {
        for (int j = 0; j < col_num; ++j) {
            x2d[i][j] = x1d[row_num * i + j];
        }
    }
    return x2d;
}

static matrixv valary1d_to_vtr2d_trans1d(valarray<double> &x1d) {
    size_t col_num = x1d.size();
    matrixv x2d(col_num, vecd(1, 0));
    for (int i = 0; i < col_num; ++i) {
        x2d[i][0] = x1d[i];
    }
    return x2d;
}

static valarray<double> vtr1d_to_valary1d(vecd &x1d) {
    size_t col_num = x1d.size();
    valarray<double> varx1d(col_num);
    for (int i = 0; i < col_num; ++i) {
        varx1d[i] = x1d[i];
    }
    return varx1d;
}



#endif //MAIN_CPP_TYPECONVERSION_H
