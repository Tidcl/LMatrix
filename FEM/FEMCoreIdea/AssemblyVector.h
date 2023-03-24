//
// Created by admin on 2023/1/6.
//

#ifndef MAIN_CPP_ASMVTR1D_H
#define MAIN_CPP_ASMVTR1D_H

#include <iostream>
#include <string>
#include "../NMatrixDefine.h"

using namespace std;

class AsmVtr1D {
    string coefficient_function_name;
    int vector_size;
    int number_of_element;
    matrixv P;
    matrixv T;
    matrixv Tb_test;
    int number_of_local_basis_fun_test;
    int basis_type_test;
    int basis_der_x_b_test;
    int basis_der_y_b_test;
    int n_gauss;
public:
    AsmVtr1D(string &coefficient_function_name, int vector_size, int number_of_element,
             matrixv &P, matrixv &T, matrixv &Tb_test,
             int number_of_local_basis_fun_test, int basis_type_test, int basis_der_x_b_test, int n_gauss);

    AsmVtr1D(string &coefficient_function_name, int vector_size, int number_of_element,
             matrixv &P, matrixv &T, matrixv &Tb_test,
             int number_of_local_basis_fun_test, int basis_type_test,
             int basis_der_x_b_test, int basis_der_y_b_test ,int n_gauss);

    ~AsmVtr1D() = default;

    matrixv assembly_vector_1d();

    matrixv assembly_vector_2d();

    double gauss_1d_test(vecd &gauss_weights, vecd &gauss_nodes,
                         matrixv &node_vertices, int &basis_index_test);

    double gauss_2d_test(vecd &gauss_weights, matrixv &gauss_nodes,
                         matrixv &node_vertices, int &basis_index_test);

};

//template<class Tv = vecd>
//Tv operator*(const double C, const Tv &x) {
//    Tv res_vtr(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) res_vtr[i] = C * x[i];
//    return res_vtr;
//}
//
//template<class Tv = vecd>
//Tv operator+(const Tv &x, const double C) {
//    Tv result_vtr(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) result_vtr[i] = C + x[i];
//    return result_vtr;
//}

//void print_matrix_asmb(matrixv &matrix) {
//    for (int i = 0; i < matrix.size(); ++i) {
//        for (int j = 0; j < matrix[0].size(); ++j) {
//            cout << matrix[i][j] << '\t' << '\t';
//        }
//        cout << endl;
//    }
//}

#endif //MAIN_CPP_ASMVTR1D_H
