//
// Created by admin on 2023/1/5.
//

#include <cmath>
#include<algorithm>
#include "AssemblyMatrix.h"
#include "FeBasisFunction.h"
#include "CoefficientFunction.h"
#include "../Tools/SymbolsOperat.h"
#include "../Tools/PrintTools.h"


using namespace std;

// 定义组装整体刚度矩阵函数，定义形参【系数函数名字，整体刚度矩阵尺寸，单元数量，局部试探基函数个数，局部测试基函数个数
//                              P矩阵，T矩阵，Tb_test矩阵 Tb_try矩阵，
//                              试探基函数类型，试探基函数阶数，测试基函数类型，测试基函数阶数】，形参基函数与刚度矩阵A中相对(笔记)
AsmMtx1D::AsmMtx1D(string &coefficient_function_name, vector<int> &matrix_size, int number_of_element,
                   int number_of_local_basis_fun_try, int number_of_local_basis_fun_test, vector<vector<double>> &P,
                   vector<vector<double>> &T, vector<vector<double>> &Tb_test, vector<vector<double>> &Tb_try,
                   int basis_type_try, int basis_type_test, int basis_der_x_try, int basis_der_x_test, int n_gauss) {
    this->coefficient_function_name = coefficient_function_name;
    this->matrix_size = matrix_size;
    this->number_of_element = number_of_element;
    this->number_of_local_basis_fun_try = number_of_local_basis_fun_try;
    this->number_of_local_basis_fun_test = number_of_local_basis_fun_test;
    this->P = P;
    this->T = T;
    this->Tb_try = Tb_try;
    this->Tb_test = Tb_test;
    this->basis_type_try = basis_type_try;
    this->basis_type_test = basis_type_test;
    this->basis_der_x_try = basis_der_x_try;
    this->basis_der_x_test = basis_der_x_test;
    this->n_gauss = n_gauss;
}

AsmMtx1D::AsmMtx1D(string &coefficient_function_name, vector<int> &matrix_size, int number_of_element,
                   int number_of_local_basis_fun_try, int number_of_local_basis_fun_test, vector<vector<double>> &P,
                   vector<vector<double>> &T, vector<vector<double>> &Tb_test, vector<vector<double>> &Tb_try,
                   int basis_type_try, int basis_type_test, int basis_der_x_try, int basis_der_x_test,
                   int basis_der_y_try, int basis_der_y_test, int n_gauss) {
    this->coefficient_function_name = coefficient_function_name;
    this->matrix_size = matrix_size;
    this->number_of_element = number_of_element;
    this->number_of_local_basis_fun_try = number_of_local_basis_fun_try;
    this->number_of_local_basis_fun_test = number_of_local_basis_fun_test;
    this->P = P;
    this->T = T;
    this->Tb_try = Tb_try;
    this->Tb_test = Tb_test;
    this->basis_type_try = basis_type_try;
    this->basis_type_test = basis_type_test;
    this->basis_der_x_try = basis_der_x_try;
    this->basis_der_x_test = basis_der_x_test;
    this->basis_der_y_try = basis_der_y_try;
    this->basis_der_y_test = basis_der_y_test;
    this->n_gauss = n_gauss;
}

// 生成整体刚度矩阵A与单元刚度矩阵S的数据结构
vector<vector<double>> AsmMtx1D::assembly_matrix_1d() {
    vector<vector<double>> vertices(1, vector<double>(T.size(), 0.0));
//    print_matrix_asm(vertices);
    vector<vector<double>> A(matrix_size[0], vector<double>(matrix_size[0], 0));
    vector<vector<double>> S(number_of_local_basis_fun_test, vector<double>(number_of_local_basis_fun_try, 0));

// 遍历所有的单元
    for (int n = 0; n < number_of_element; ++n) {
//      获取每个单元的全局节点坐标，先通过T找出第n的单元列，
//      代表了第n个单元的节点编号，再通过P找出第n个单元对应节点编号的对应坐标
//        vertices = P[:,T[:,n]];
        for (int i = 0; i < P.size(); ++i) {
            for (int j = 0; j < T.size(); ++j) {
                for (int k = 0; k < n + 1; ++k) {
                    vertices[i][j] = P[i][int(T[j][k])];
                }
            }
        }
//        print_matrix_asm(vertices);
//      调用方法，传入高斯点数目，获取区间[-1,1]上的高斯点位置与权数
        vector<double> gauss_nodes_reference;
        vector<double> gauss_weights_reference;
        vector<double> gauss_nodes_local;
        vector<double> gauss_weights_local;

        if (n_gauss == 3) {
            gauss_nodes_reference = {-0.77459667, 0.0, 0.77459667};
            gauss_weights_reference = {0.55555556, 0.888888889, 0.55555556};
        } else if (n_gauss == 4) {
            gauss_nodes_reference = {-0.86113631, -0.33998104, 0.33998104, 0.86113631};
            gauss_weights_reference = {0.34785485, 0.65214515, 0.65214515, 0.34785485};
        } else if (n_gauss == 5) {
            gauss_nodes_reference = {-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985};
            gauss_weights_reference = {0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689};
        } else {
            cout << "Error: I haven't written it yet!" << endl;
        }
//      修改标准高斯积分的坐标值和点值，高斯积分在标准区间的权值和位置点如上，
//      将其修改为单元积分区间的权值和位置点，和载荷向量处是重复的
        double upper_limit = *max_element(vertices[0].begin(), vertices[0].end());
        double lower_limit = *min_element(vertices[0].begin(), vertices[0].end());
        gauss_nodes_local = (upper_limit - lower_limit) / 2 * gauss_nodes_reference + (upper_limit + lower_limit) / 2;
        gauss_weights_local = (upper_limit - lower_limit) / 2 * gauss_weights_reference;

//      全局基函数的下角标i，j，替换为局部基函数alpha，beta。(i，j也表示了两个基函数，测试与试探)
        for (int alpha = 0; alpha < number_of_local_basis_fun_try; ++alpha) {
            for (int beta = 0; beta < number_of_local_basis_fun_test; ++beta) {
//              调用了在刚度矩阵中计算积分值r的具体函数，传实参【】，有的在本函数体给出，vertices等，有的通过形参往上层函数调，有点容易乱
                double int_value = gauss_1d_try_test(gauss_weights_local, gauss_nodes_local,
                                                     vertices, alpha, beta);
                S[beta][alpha] = int_value;  // 单元刚度
                int i = int(Tb_test[beta][n]);  // b的列数是第n个单元的编号，因此：第n个单元第beta个节点的全局编号
                int j = int(Tb_try[alpha][n]);  // 第n个单元第alpha个节点的全局编号
                A[i][j] += int_value;  // 组装整体刚度矩阵
            }
        }
    }
//    print_matrix_asm(A);
    return A;
}

//  定义了在刚度矩阵中计算积分值r的具体函数，定义形参【系数函数名称，高斯点权数，高斯点位置，单元节点坐标，
//                                           试探基函数类型，判断试探与测试基函数的索引，试探基函数阶数
//                                           测试基函数类型，判断试探与测试基函数的索引，测试基函数阶数】
double AsmMtx1D::gauss_1d_try_test(vector<double> &gauss_weights, vector<double> &gauss_nodes,
                                   vector<vector<double>> &node_vertices, int &basis_index_try, int &basis_index_test) {
//  说明高斯点数目
    int gpn = n_gauss;
//  正儿八经计算刚度矩阵中的积分值int_value
    double int_value = 0.0;

//  通过CoefficientFun类生成实例的时候传入选择标记，从而获取系数函数
    CoefficientFun coefficientFun('c');
    double cof_func = 0;
//  求值时，调用生成局部基函数的函数，传实参【】，最后作为形参往上层函数调，有点容易乱，
//  在这里试探/测试基函数用的同一个'fe_basis_fun_1d'
    for (int i = 0; i < gpn; ++i) {
        BasisFun1D tryBasisFun1D(gauss_nodes[i], node_vertices, basis_type_try, basis_index_try, basis_der_x_try);
        BasisFun1D testBasisFun1D(gauss_nodes[i], node_vertices, basis_type_test, basis_index_test, basis_der_x_test);
        int_value += gauss_weights[i] * coefficientFun.coefficient_function(gauss_nodes[i], cof_func)
                     * tryBasisFun1D.fe_basis_fun_1d() * testBasisFun1D.fe_basis_fun_1d();
    }
    return int_value;
}


vector<vector<double>> AsmMtx1D::assembly_matrix_2d() {
    vector<vector<double>> vertices(2, vector<double>(T.size(), 0.0));
//    print_matrix_asm(vertices);
    vector<vector<double>> A(matrix_size[0], vector<double>(matrix_size[0], 0));
    vector<vector<double>> S(number_of_local_basis_fun_test, vector<double>(number_of_local_basis_fun_try, 0));

//    vertices = P[:,T[:,n]];
    for (int n = 0; n < number_of_element; ++n) {
        for (int i = 0; i < P.size(); ++i) {
            for (int j = 0; j < T.size(); ++j) {
                for (int k = 0; k < n + 1; ++k) {
                    vertices[i][j] = P[i][int(T[j][k])];
                }
            }
        }

        vector<vector<double>> gauss_nodes_reference = {{0.5, 0},
                                                        {0.5, 0.5},
                                                        {0,   0.5}};
        vector<double> gauss_weights_reference = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};

        double x1 = vertices[0][0];
        double y1 = vertices[1][0];
        double x2 = vertices[0][1];
        double y2 = vertices[1][1];
        double x3 = vertices[0][2];
        double y3 = vertices[1][2];
        double jacobi = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

        vector<double> gauss_weights_local = jacobi * gauss_weights_reference;
        vector<vector<double>> gauss_nodes_local(3, vector<double>(2, 0.0));


        for (size_t m = 0; m < gauss_nodes_local.size(); ++m) {
            gauss_nodes_local[m][0] =
                    x1 + (x2 - x1) * gauss_nodes_reference[m][0] + (x3 - x1) * gauss_nodes_reference[m][1];
            gauss_nodes_local[m][1] =
                    y1 + (y2 - y1) * gauss_nodes_reference[m][0] + (y3 - y1) * gauss_nodes_reference[m][1];
        }

        for (int alpha = 0; alpha < number_of_local_basis_fun_try; ++alpha) {
            for (int beta = 0; beta < number_of_local_basis_fun_test; ++beta) {
//              调用了在刚度矩阵中计算积分值r的具体函数，传实参【】，有的在本函数体给出，vertices等，有的通过形参往上层函数调，有点容易乱
                double int_value = gauss_2d_try_test(gauss_weights_local, gauss_nodes_local,
                                                     vertices, alpha, beta);
                S[beta][alpha] = int_value;  // 单元刚度
                int i = int(Tb_test[beta][n]);  // b的列数是第n个单元的编号，因此：第n个单元第beta个节点的全局编号
                int j = int(Tb_try[alpha][n]);  // 第n个单元第alpha个节点的全局编号
                A[i][j] += int_value;  // 组装整体刚度矩阵
            }
        }
    }
    return A;
}

double AsmMtx1D::gauss_2d_try_test(vector<double> &gauss_weights, vector<vector<double>> &gauss_nodes,
                                   vector<vector<double>> &node_vertices, int &basis_index_try, int &basis_index_test) {

    int gpn = n_gauss;
    double int_value = 0.0;
    CoefficientFun coefficientFun('c');
    double cof_func = 0;
    double temp_cof_y = 0;

    for (int i = 0; i < gpn; ++i) {
        BasisFun1D tryBasisFun2D(gauss_nodes[i][0], gauss_nodes[i][1], node_vertices,
                                 basis_type_try, basis_index_try, basis_der_x_try, basis_der_y_try);
        BasisFun1D testBasisFun2D(gauss_nodes[i][0], gauss_nodes[i][1], node_vertices,
                                  basis_type_test, basis_index_test, basis_der_x_test, basis_der_y_test);

        int_value += gauss_weights[i] * coefficientFun.coefficient_function_2d(gauss_nodes[i][0], temp_cof_y, cof_func)
                     * coefficientFun.coefficient_function_2d(gauss_nodes[i][1], temp_cof_y, cof_func)
                     * tryBasisFun2D.fe_local_basis_fun_2d() * testBasisFun2D.fe_local_basis_fun_2d();
    }
    return int_value;

}

//int main() {
//    //test 2 pass
//    string str_test = "niulirui";
//    vector<int> mtx_size = {5, 5};
//    int num_of_ele = 4;
//    int number_of_local_basis_fun_try = 2;
//    int number_of_local_basis_fun_test = 2;
//    vector<vector<double>> P = {{0, 0.25, 0.5, 0.75, 1}};
//    vector<vector<double>> T = {{0, 1, 2, 3},
//                                {1, 2, 3, 4}};
//    vector<vector<double>> Tb_test = {{0, 1, 2, 3},
//                                      {1, 2, 3, 4}};
//    vector<vector<double>> Tb_try = {{0, 1, 2, 3},
//                                     {1, 2, 3, 4}};
//    int basis_type_try = 101;
//    int basis_type_test = 101;
//    int basis_der_x_try = 1;
//    int basis_der_x_test = 1;
//    int n_gauss = 4;
//
//    AsmMtx1D test_fun_2(str_test, mtx_size, num_of_ele, number_of_local_basis_fun_try, number_of_local_basis_fun_test,
//                        P, T, Tb_test, Tb_try, basis_type_try, basis_type_test, basis_der_x_try, basis_der_x_test,
//                        n_gauss);
//
//    test_fun_2.assembly_matrix_1d();
//
//    cout << "hahaha" << endl;
//}

//int main() {
//    //test 2 pass
//    string str_test = "c";
//    vector<int> mtx_size = {9, 9};
//    int num_of_ele = 8;
//    int number_of_local_basis_fun_try = 3;
//    int number_of_local_basis_fun_test = 3;
//    vector<vector<double>> P = {{-1.0, -1.0, -1.0, 0.0,  0.0, 0.0, 1.0,  1.0, 1.0},
//                                {-1.0, 0.0,  1.0,  -1.0, 0.0, 1.0, -1.0, 0.0, 1.0}};
//    vector<vector<double>> T = {{0, 1, 1, 2, 3, 4, 4, 5},
//                                {3, 3, 4, 4, 6, 6, 7, 7},
//                                {1, 4, 2, 5, 4, 7, 5, 8}};
//    vector<vector<double>> Tb_test = {{0, 1, 1, 2, 3, 4, 4, 5},
//                                      {3, 3, 4, 4, 6, 6, 7, 7},
//                                      {1, 4, 2, 5, 4, 7, 5, 8}};
//    vector<vector<double>> Tb_try = {{0, 1, 1, 2, 3, 4, 4, 5},
//                                     {3, 3, 4, 4, 6, 6, 7, 7},
//                                     {1, 4, 2, 5, 4, 7, 5, 8}};
//    int basis_type_try = 201;
//    int basis_type_test = 201;
//
//    int basis_der_x_try = 1;
//    int basis_der_y_try = 0;
//    int basis_der_x_test = 1;
//    int basis_der_y_test = 0;
//
//    int n_gauss = 3;
//
//    AsmMtx1D test_fun_2(str_test, mtx_size, num_of_ele,
//                        number_of_local_basis_fun_try, number_of_local_basis_fun_test,
//                        P, T, Tb_test, Tb_try,
//                        basis_type_try, basis_type_test,
//                        basis_der_x_try, basis_der_x_test,
//                        basis_der_y_try, basis_der_y_test, n_gauss);
//
//    test_fun_2.assembly_matrix_2d();
//
//    cout << "hahaha" << endl;
//}