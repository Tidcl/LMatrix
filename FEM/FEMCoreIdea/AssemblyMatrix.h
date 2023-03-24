//
// Created by admin on 2023/1/5.
//

#ifndef MAIN_CPP_ASMMTX1D_H
#define MAIN_CPP_ASMMTX1D_H

#include <iostream>
#include <string>
#include "../NMatrixDefine.h"

using namespace std;
class AsmMtx1D {
private:
    string coefficient_function_name;
    veci matrix_size;
    int number_of_element;
    int number_of_local_basis_fun_try;
    int number_of_local_basis_fun_test;
    matrixv P;
    matrixv T;
    matrixv Tb_try;
    matrixv Tb_test;
    int basis_type_try;
    int basis_type_test;
    int basis_der_x_try;
    int basis_der_x_test;
    int basis_der_y_try{};
    int basis_der_y_test{};
    int n_gauss;
public:
    AsmMtx1D(string &coefficient_function_name, veci &matrix_size, int number_of_element,
             int number_of_local_basis_fun_try, int number_of_local_basis_fun_test,
             matrixv &P, matrixv &T,
             matrixv &Tb_test, matrixv &Tb_try,
             int basis_type_try, int basis_type_test,
             int basis_der_x_try, int basis_der_x_test, int n_gauss);

    AsmMtx1D(string &coefficient_function_name, veci &matrix_size, int number_of_element,
             int number_of_local_basis_fun_try, int number_of_local_basis_fun_test,
             matrixv &P, matrixv &T,
             matrixv &Tb_test, matrixv &Tb_try,
             int basis_type_try, int basis_type_test,
             int basis_der_x_try, int basis_der_x_test,
             int basis_der_y_try, int basis_der_y_test, int n_gauss);

    ~AsmMtx1D() = default;

    matrixv assembly_matrix_1d();

    matrixv assembly_matrix_2d();

    double gauss_1d_try_test(vecd &gauss_weights, vecd &gauss_nodes,
                             matrixv &node_vertices, int &basis_index_try, int &basis_index_test);

    double gauss_2d_try_test(vecd &gauss_weights, matrixv &gauss_nodes,
                             matrixv &node_vertices, int &basis_index_try, int &basis_index_test);
};

//重载运算符的时候，需要注意运算符左右两侧的数据类型的顺序
//template<class Tv = vecd>
//Tv operator*(const double C, const Tv &x) {
//    Tv res_mtx(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) res_mtx[i] = C * x[i];
//    return res_mtx;
//}
//
//
//template<class Tv = vecd>
//Tv operator+(const Tv &x, const double C) {
//    Tv result_mtx(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) result_mtx[i] = C + x[i];
//    return result_mtx;
//}

//void print_matrix_asm(matrixv &matrix) {
//    for (int i = 0; i < matrix.size(); ++i) {
//        for (int j = 0; j < matrix[0].size(); ++j) {
//            cout << matrix[i][j] << '\t' << '\t';
//        }
//        cout << endl;
//    }
//}

#endif //MAIN_CPP_ASMMTX1D_H
