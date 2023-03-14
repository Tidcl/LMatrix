//
// Created by admin on 2023/1/5.
//

#ifndef MAIN_CPP_ASMMTX1D_H
#define MAIN_CPP_ASMMTX1D_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class AsmMtx1D {
private:
    string coefficient_function_name;
    vector<int> matrix_size;
    int number_of_element;
    int number_of_local_basis_fun_try;
    int number_of_local_basis_fun_test;
    vector<vector<double>> P;
    vector<vector<double>> T;
    vector<vector<double>> Tb_try;
    vector<vector<double>> Tb_test;
    int basis_type_try;
    int basis_type_test;
    int basis_der_x_try;
    int basis_der_x_test;
    int basis_der_y_try{};
    int basis_der_y_test{};
    int n_gauss;
public:
    AsmMtx1D(string &coefficient_function_name, vector<int> &matrix_size, int number_of_element,
             int number_of_local_basis_fun_try, int number_of_local_basis_fun_test,
             vector<vector<double>> &P, vector<vector<double>> &T,
             vector<vector<double>> &Tb_test, vector<vector<double>> &Tb_try,
             int basis_type_try, int basis_type_test,
             int basis_der_x_try, int basis_der_x_test, int n_gauss);

    AsmMtx1D(string &coefficient_function_name, vector<int> &matrix_size, int number_of_element,
             int number_of_local_basis_fun_try, int number_of_local_basis_fun_test,
             vector<vector<double>> &P, vector<vector<double>> &T,
             vector<vector<double>> &Tb_test, vector<vector<double>> &Tb_try,
             int basis_type_try, int basis_type_test,
             int basis_der_x_try, int basis_der_x_test,
             int basis_der_y_try, int basis_der_y_test, int n_gauss);

    ~AsmMtx1D() = default;

    vector<vector<double>> assembly_matrix_1d();

    vector<vector<double>> assembly_matrix_2d();

    double gauss_1d_try_test(vector<double> &gauss_weights, vector<double> &gauss_nodes,
                             vector<vector<double>> &node_vertices, int &basis_index_try, int &basis_index_test);

    double gauss_2d_try_test(vector<double> &gauss_weights, vector<vector<double>> &gauss_nodes,
                             vector<vector<double>> &node_vertices, int &basis_index_try, int &basis_index_test);
};

//重载运算符的时候，需要注意运算符左右两侧的数据类型的顺序
//template<class Tv = vector<double>>
//Tv operator*(const double C, const Tv &x) {
//    Tv res_mtx(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) res_mtx[i] = C * x[i];
//    return res_mtx;
//}
//
//
//template<class Tv = vector<double>>
//Tv operator+(const Tv &x, const double C) {
//    Tv result_mtx(x.size(), 0);
//    for (int i = 0; i < x.size(); ++i) result_mtx[i] = C + x[i];
//    return result_mtx;
//}

//void print_matrix_asm(vector<vector<double>> &matrix) {
//    for (int i = 0; i < matrix.size(); ++i) {
//        for (int j = 0; j < matrix[0].size(); ++j) {
//            cout << matrix[i][j] << '\t' << '\t';
//        }
//        cout << endl;
//    }
//}

#endif //MAIN_CPP_ASMMTX1D_H
