//
// Created by Administrater on 2023/1/4.
//

#ifndef RUNSCRIPY_PY_GENERATEMATRIXPBTB_H
#define RUNSCRIPY_PY_GENERATEMATRIXPBTB_H

#include <iostream>
#include "../NMatrixDefine.h"

using namespace std;

class GenPbTb {
private:
    double left, right, bottom{},top{},h{};
    int basis_type_try;
    int basis_type_test;
    vecd h2;
public:
    GenPbTb(double left, double right, double h, int basis_type_try, int basis_type_test);
    GenPbTb(double left, double right,double bottom,double top, vecd&h2, int basis_type_try, int basis_type_test);

    ~GenPbTb() = default;

    void gen_pb_tb_1d(matrixv &P, matrixv &T, double &h_basis) const;
    void gen_pb_tb_2d(matrixv &P,matrixv &T);
};

//void print_vector(vecd &array) {
//    for (int i = 0; i < array.size(); ++i) {
//        cout << array[i] << '\t' << '\t';
//    }
//    cout << endl;
//}

//void print_matrix(matrixv &matrix) {
//    for (int i = 0; i < matrix.size(); ++i) {
//        for (int j = 0; j < matrix[0].size(); ++j) {
//            cout << matrix[i][j] << '\t' << '\t';
//        }
//        cout << endl;
//    }
//}


#endif //RUNSCRIPY_PY_GENERATEMATRIXPBTB_H
