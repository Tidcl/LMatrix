//
// Created by admin on 2023/1/6.
//
#include <algorithm>
#include "AssemblyVector.h"
#include "CoefficientFunction.h"
#include "FeBasisFunction.h"
#include "../Tools/SymbolsOperat.h"
#include "../Tools/PrintTools.h"


// ***此部分与上文装配整体刚度矩阵A类似，因此注释简略***————区别在于 没有试探基函数-None-try
// 定义组装载荷向量函数，定义形参【系数函数名字，载荷向量尺寸，单元数量，P矩阵，T矩阵，Tb_test矩阵，
//                           局部测试基函数个数，测试基函数类型，测试基函数阶段】，形参基函数与载荷向量b中相对(笔记)
AsmVtr1D::AsmVtr1D(std::string &coefficient_function_name, int vector_size, int number_of_element,
                   vector<vector<double>> &P, vector<vector<double>> &T, vector<vector<double>> &Tb_test,
                   int number_of_local_basis_fun_test, int basis_type_test, int basis_der_x_b_test, int n_gauss) {
    this->coefficient_function_name = coefficient_function_name;
    this->vector_size = vector_size;
    this->number_of_element = number_of_element;
    this->P = P;
    this->T = T;
    this->Tb_test = Tb_test;
    this->number_of_local_basis_fun_test = number_of_local_basis_fun_test;
    this->basis_type_test = basis_type_test;
    this->basis_der_x_b_test = basis_der_x_b_test;
    this->n_gauss = n_gauss;
}

AsmVtr1D::AsmVtr1D(std::string &coefficient_function_name, int vector_size, int number_of_element,
                   vector<vector<double>> &P, vector<vector<double>> &T, vector<vector<double>> &Tb_test,
                   int number_of_local_basis_fun_test, int basis_type_test,
                   int basis_der_x_b_test, int basis_der_y_b_test, int n_gauss) {
    this->coefficient_function_name = coefficient_function_name;
    this->vector_size = vector_size;
    this->number_of_element = number_of_element;
    this->P = P;
    this->T = T;
    this->Tb_test = Tb_test;
    this->number_of_local_basis_fun_test = number_of_local_basis_fun_test;
    this->basis_type_test = basis_type_test;
    this->basis_der_x_b_test = basis_der_x_b_test;
    this->basis_der_y_b_test = basis_der_y_b_test;
    this->n_gauss = n_gauss;
}

vector<vector<double>> AsmVtr1D::assembly_vector_1d() {
//  生成载荷向量b的数据结构
    vector<vector<double>> b(vector_size, vector<double>(1, 0));

    vector<vector<double>> vertices(1, vector<double>(T.size(), 0));
    for (int n = 0; n < number_of_element; ++n) {
//      vertices = P[:,T[:,n]];
        for (int i = 0; i < P.size(); ++i) {
            for (int j = 0; j < T.size(); ++j) {
                for (int k = 0; k < n + 1; ++k) {
                    vertices[i][j] = P[i][int(T[j][k])];
                }
            }
        }
//        print_matrix_asmb(vertices);

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

        double upper_limit = *max_element(vertices[0].begin(), vertices[0].end());
        double lower_limit = *min_element(vertices[0].begin(), vertices[0].end());
        gauss_nodes_local = (upper_limit - lower_limit) / 2 * gauss_nodes_reference + (upper_limit + lower_limit) / 2;
        gauss_weights_local = (upper_limit - lower_limit) / 2 * gauss_weights_reference;

//      只有一个基函数对应下角标i，对应为局部基函数，beta。(只有一个测试基函数)
        for (int beta = 0; beta < number_of_local_basis_fun_test; ++beta) {
            double temp = gauss_1d_test(gauss_weights_local, gauss_nodes_local, vertices, beta);
            int i = int(Tb_test[beta][n]);  // 同上解释，第n个单元第beta个节点的全局编号
            b[i][0] += temp;  // 组装载荷向量
        }
    }
//    myPrint.print_matrix(b);
//    print_matrix_asmb(b);
    return b;
}


// ***此部分与上文在刚度矩阵中计算积分值r类似，因此注释简略***————区别在于没有试探基函数
// 定义了在载荷向量中计算积分值r的具体函数，定义形参【系数函数名称，高斯点权数，高斯点位置，单元节点坐标，
//                                          测试基函数类型，判断试探与测试基函数的索引，测试基函数阶数】
double AsmVtr1D::gauss_1d_test(vector<double> &gauss_weights, vector<double> &gauss_nodes,
                               vector<vector<double>> &node_vertices, int &basis_index_test) {

    int gpn = n_gauss;
    double int_result = 0.0;

    CoefficientFun coefficientFun('f');
    double cof_func = 0;
//  正儿八经计算载荷向量中的积分值int_result，解释同上
    for (int i = 0; i < gpn; ++i) {
        BasisFun1D basisFun1D(gauss_nodes[i], node_vertices, basis_type_test, basis_index_test, basis_der_x_b_test);
        int_result += gauss_weights[i] * coefficientFun.coefficient_function(gauss_nodes[i], cof_func)
                      * basisFun1D.fe_basis_fun_1d();
    }
    return int_result;
}

vector<vector<double>> AsmVtr1D::assembly_vector_2d() {
    //  生成载荷向量b的数据结构
    vector<vector<double>> b(vector_size, vector<double>(1, 0));

    vector<vector<double>> vertices(2, vector<double>(T.size(), 0));
    for (int n = 0; n < number_of_element; ++n) {
//      vertices = P[:,T[:,n]];
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

//      只有一个基函数对应下角标i，对应为局部基函数，beta。(只有一个测试基函数)
        for (int beta = 0; beta < number_of_local_basis_fun_test; ++beta) {
            double temp = gauss_2d_test(gauss_weights_local, gauss_nodes_local, vertices, beta);
            int i = int(Tb_test[beta][n]);  // 同上解释，第n个单元第beta个节点的全局编号
            b[i][0] += temp;  // 组装载荷向量
        }
    }
//    int sds = 8;
//    MyPrint myPrint(sds);
//    myPrint.print_matrix(b);
    return b;
//    [[ 0.13015926]
//    [ 0.01528012]
//    [-0.7346996 ]
//    [ 0.2613233 ]
//    [-0.26371372]
//    [ 0.40286974]
//    [-0.16864344]
//    [ 0.49091957]
//    [ 0.93368522]]
}

double AsmVtr1D::gauss_2d_test(vector<double> &gauss_weights, vector<vector<double>> &gauss_nodes,
                               vector<vector<double>> &node_vertices, int &basis_index_test) {

    int gpn = n_gauss;
    double int_result = 0.0;

    CoefficientFun coefficientFun('f');
    double cof_func = 0;

    for (int i = 0; i < gpn; ++i) {
        BasisFun1D basisFun2D(gauss_nodes[i][0], gauss_nodes[i][1], node_vertices, basis_type_test,
                              basis_index_test, basis_der_x_b_test, basis_der_y_b_test);
        int_result += gauss_weights[i] *
                      coefficientFun.coefficient_function_2d(gauss_nodes[i][0], gauss_nodes[i][1], cof_func)
                      * basisFun2D.fe_local_basis_fun_2d();
    }
    return int_result;
}
//int main() {
//    //test 3 pass
//    string str_test = "niulirui";
//    int vtr_size = 5;
//    int num_of_ele = 4;
//    vector<vector<double>> P = {{0, 0.25, 0.5, 0.75, 1}};
//    vector<vector<double>> T = {{0, 1, 2, 3},
//                                {1, 2, 3, 4}};
//    vector<vector<double>> Tb_test = {{0, 1, 2, 3},
//                                      {1, 2, 3, 4}};
//    int number_of_local_basis_fun_test = 2;
//    int basis_type_test = 101;
//    int basis_der_x_test = 0;
//    int n_gauss = 4;
//
//
//    AsmVtr1D test_fun_3(str_test, vtr_size, num_of_ele, P, T, Tb_test, number_of_local_basis_fun_test, basis_type_test,
//                        basis_der_x_test, n_gauss);
//
//    test_fun_3.assembly_vector_1d();
//
//    cout << "hahaha" << endl;
//}

//int main() {
//    //test 3 pass
//    string str_test = "f";
//    int vtr_size = 9;
//    int num_of_ele = 8;
//    vector<vector<double>> P = {{-1.0, -1.0, -1.0, 0.0,  0.0, 0.0, 1.0,  1.0, 1.0},
//                                {-1.0, 0.0,  1.0,  -1.0, 0.0, 1.0, -1.0, 0.0, 1.0}};
//    vector<vector<double>> T = {{0, 1, 1, 2, 3, 4, 4, 5},
//                                {3, 3, 4, 4, 6, 6, 7, 7},
//                                {1, 4, 2, 5, 4, 7, 5, 8}};
//    vector<vector<double>> Tb_test = {{0, 1, 1, 2, 3, 4, 4, 5},
//                                      {3, 3, 4, 4, 6, 6, 7, 7},
//                                      {1, 4, 2, 5, 4, 7, 5, 8}};
//    int number_of_local_basis_fun_test = 3;
//    int basis_type_test = 201;
//    int basis_der_x_test = 0;
//    int basis_der_y_test = 0;
//    int n_gauss = 3;
//
//
//    AsmVtr1D test_fun_3(str_test, vtr_size, num_of_ele, P, T, Tb_test, number_of_local_basis_fun_test, basis_type_test,
//                        basis_der_x_test,basis_der_y_test, n_gauss);
//
//    test_fun_3.assembly_vector_2d();
//
//    cout << "hahaha" << endl;
//}