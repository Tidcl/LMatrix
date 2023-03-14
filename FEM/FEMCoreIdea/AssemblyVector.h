//
// Created by admin on 2023/1/6.
//

#ifndef MAIN_CPP_ASMVTR1D_H
#define MAIN_CPP_ASMVTR1D_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class AsmVtr1D {
    string coefficient_function_name;
    int vector_size;
    int number_of_element;
    vector<vector<double>> P;
    vector<vector<double>> T;
    vector<vector<double>> Tb_test;
    int number_of_local_basis_fun_test;
    int basis_type_test;
    int basis_der_x_b_test;
    int basis_der_y_b_test;
    int n_gauss;
public:
    AsmVtr1D(string &coefficient_function_name, int vector_size, int number_of_element,
             vector<vector<double>> &P, vector<vector<double>> &T, vector<vector<double>> &Tb_test,
             int number_of_local_basis_fun_test, int basis_type_test, int basis_der_x_b_test, int n_gauss);

    AsmVtr1D(string &coefficient_function_name, int vector_size, int number_of_element,
             vector<vector<double>> &P, vector<vector<double>> &T, vector<vector<double>> &Tb_test,
             int number_of_local_basis_fun_test, int basis_type_test,
             int basis_der_x_b_test, int basis_der_y_b_test ,int n_gauss);

    ~AsmVtr1D() = default;

    vector<vector<double>> assembly_vector_1d();

    vector<vector<double>> assembly_vector_2d();

    double gauss_1d_test(vector<double> &gauss_weights, vector<double> &gauss_nodes,
                         vector<vector<double>> &node_vertices, int &basis_index_test);

    double gauss_2d_test(vector<double> &gauss_weights, vector<vector<double>> &gauss_nodes,
                         vector<vector<double>> &node_vertices, int &basis_index_test);

};

//template<class Tv = vector<double>>
//Tv operator*(const double C, const Tv &x) {
//    Tv res_vtr(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) res_vtr[i] = C * x[i];
//    return res_vtr;
//}
//
//template<class Tv = vector<double>>
//Tv operator+(const Tv &x, const double C) {
//    Tv result_vtr(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) result_vtr[i] = C + x[i];
//    return result_vtr;
//}

//void print_matrix_asmb(vector<vector<double>> &matrix) {
//    for (int i = 0; i < matrix.size(); ++i) {
//        for (int j = 0; j < matrix[0].size(); ++j) {
//            cout << matrix[i][j] << '\t' << '\t';
//        }
//        cout << endl;
//    }
//}

#endif //MAIN_CPP_ASMVTR1D_H
