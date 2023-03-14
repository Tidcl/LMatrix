//
// Created by Administrater on 2023/1/7.
//

#include <sys/time.h>
#include "../Tools/PrintTools.h"
#include "../Tools/TypeConversion.h"
#include "../Tools/SymbolsOperat.h"
#include "../SolverMethod/IterativeMethod.h"
#include "../SolverMethod/DirectMethod.h"
#include "../SolverMethod/SolveEigenValue.h"
#include "../ErrorCompare/ErrorCompare.h"
#include "FEMSolver.h"
#include "GenerateMatrixPbTb.h"
#include "AssemblyMatrix.h"
#include "AssemblyVector.h"
#include "BoundaryInformation.h"


// try试探为真物理，test测试为人造的。全局下，1个单元，1组基函数，局部下，对应在不同节点上有不同的基函数分量(个数)，
// 定义了1维求解函数的框架类，定义类初始变量【左侧点坐标（1个），右侧点坐标（1个），尺寸，试探基(形)函数类型，测试基(形)函数类型，积分点数量】
FeSolverPoisson1D::FeSolverPoisson1D(int left, int right, vector<double> &h,
                                     int basis_type_try, int basis_type_test,
                                     int number_of_integral) {
    this->left = left;
    this->right = right;
    this->h = h;
    this->basis_type_try = basis_type_try;
    this->basis_type_test = basis_type_test;
    this->number_of_integral = number_of_integral;
}

FeSolverPoisson1D::FeSolverPoisson1D(int left, int right, int bottom, int top, vector<double> &h, int basis_type_try,
                                     int basis_type_test, int number_of_integral) {
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->h = h;
    this->basis_type_try = basis_type_try;
    this->basis_type_test = basis_type_test;
    this->number_of_integral = number_of_integral;

}

// 框架类的成员函数，控制输出，【输出精度，输出总刚A，输出特征值，完整特征值/向量输出，有限元求解结果，误差列表】
void FeSolverPoisson1D::fe_solver_poisson_1d(int &accuracy, bool print_A = false, bool print_eigen = false,
                                             bool entire = false, bool fem_result = false, bool error_ary = false,
                                             bool fast_method = true) {
//  实例化工具中的打印对象
    MyPrint myPrint(accuracy);
//  生成网格，网格数 == 单元个数，N_mesh = 4
    int N_mesh = int((right - left) / h[0]);

    cout << "--------Global--------" << endl;
    cout << "Number of mesh = " << N_mesh << endl;


// 调用生成网格信息矩阵函数
// P：所有网格节点的坐标(1个网格基本是2个节点，可以增加，存放在P)
// T：所有网格全局编号，列数代表第i个网格，某列即代表了某个网格的顺序节点编号
// P，T是Pb，Tb的特殊情况
    vector<vector<double>> P;
    vector<vector<double>> T;
    double h_basis;

    GenPbTb genPbTb(left, right, h[0], basis_type_try, basis_type_test);
    genPbTb.gen_pb_tb_1d(P, T, h_basis);

// 获取单元个数，T是二维Vector，size()可以拿到[行，列]，列数即为单元个数，储存浮点，转换为整型作为后续的下标
    int number_of_element = int(T[0].size());
    int N_basis = int(P[0].size()) - 1;
    cout << "Number of elements = " << number_of_element << endl;
    cout << "Number of nodes = " << N_basis+1 << endl;
    cout << "--------------------------------------------------------" << endl;

//  根据试探基函数的类型，选择Pb和Tb，这里是线性，标号101，所以局部试探基函数的个数是2，
//  可以理解为：每个1维单元上有2节点，每个节点上1个基函数
    int number_of_local_basis_fun_try;
    vector<vector<double>> Tb_try;
    if (basis_type_try == 101) {
//        Pb_try = P;
        Tb_try = T;
        number_of_local_basis_fun_try = 2;
    } else if (basis_type_try == 102) {
//        Pb_try = P;
        Tb_try = T;
        number_of_local_basis_fun_try = 3;
    } else {
        cout << "Are you drinking a fake, Just 101 and 102" << endl;
        system("pause");
        exit(1);
    }


    int number_of_local_basis_fun_test;
    vector<vector<double>> Tb_test;
    if (basis_type_test == 101) {
//        Pb_test = P;
        Tb_test = T;
        number_of_local_basis_fun_test = 2;
    } else if (basis_type_test == 102) {
//        Pb_test = P;
        Tb_test = T;
        number_of_local_basis_fun_test = 3;
    } else {
        cout << "Are you drinking a fake, Just 101 and 102" << endl;
        system("pause");
        exit(1);
    }


//  指定整体刚度矩阵的行数与列数，这里是线性，一共有网格数+1的节点
    vector<int> matrix_size(N_basis + 1, N_basis + 1);
//  指定载荷向量的列数，这里是线性，一共有网格数+1的节点，与矩阵对应
    int vector_size = N_basis + 1;

//  整体刚度矩阵中 试探与测试基函数的阶数
    int basis_der_x_try = 1;
    int basis_der_x_test = 1;

//  载荷向量中，测试函数的阶数
    int basis_der_x_b_test = 0;

//  高斯积分点数
    int number_of_gauss_point = int(number_of_integral);

    string str_test = "NiuLiuRui";

//  调用组装整体刚度矩阵函数，传实参【系数矩阵拉姆达表达式，矩阵尺寸元组，单元个数整型，试探局部基函数个数整型，测试局部基函数个数整型，
//                             P矩阵二维Vector，T矩阵二维Vector，Tb_test矩阵二维Vector，Tb_try矩阵二维Vector，
//                             刚度矩阵试探基函数类型整型，刚度矩阵试探基函数阶数整型，
//                             刚度矩阵测试基函数类型整型，刚度矩阵测试基函数阶数整型】  !!!二次不变
    AsmMtx1D asmMtx1D(str_test, matrix_size, number_of_element, number_of_local_basis_fun_try,
                      number_of_local_basis_fun_test, P, T, Tb_test, Tb_try,
                      basis_type_try, basis_type_test, basis_der_x_try, basis_der_x_test,
                      number_of_gauss_point);

    vector<vector<double>> A = asmMtx1D.assembly_matrix_1d();
    if (print_A) myPrint.print_matrix(A);
//    myPrint.print_matrix(A);

//  调用组装载荷向量矩阵函数，传实参【系数矩阵拉姆达表达式，向量尺寸整型，单元个数整型，
//                              P矩阵二维Vector，T矩阵二维Vector，Tb_test矩阵二维Vector，
//                             测试局部基函数个数整型，载荷向量测试基函数类型整型，载荷向量测试基函数阶数整型】
    AsmVtr1D asmVtr1D(str_test, vector_size, number_of_element, P, T, Tb_test,
                      number_of_local_basis_fun_test, basis_type_test, basis_der_x_b_test, number_of_gauss_point);

    vector<vector<double>> b = asmVtr1D.assembly_vector_1d();
//    myPrint.print_matrix(b);

// 调用生成边界节点信息矩阵函数，传实参【单元个数整型】
    vector<vector<double>> boundary_nodes(0, vector<double>(0, 0));
    BoundaryInformation boundaryInformation(N_basis, boundary_nodes);
    boundaryInformation.generate_boundary_nodes();
// A, b, P放在函数方便引用改变形参
// 调用处理边界条件函数，传实参【未修改的整体刚度矩阵，未修改的载荷向量矩阵，边界节点信息矩阵】
    boundaryInformation.treat_boundary_condition(A, b, P, 1);
//    myPrint.print_matrix(A);
//    myPrint.print_matrix(b);

//  计算方程组的时候，A期待2维Vector，b期待1维valarray，而到目前为止，A，b都是二维的Vector，因此需要转换
    valarray<double> valary_b = vtr2d_to_valary1d(b);
//    myPrint.print_vector(valary_b);

//    Tdma tdmaSolve(valary_b.size(), A, valary_b);
//    valarray<double> valary_sol = tdmaSolve.solver();

//    GaussEpp gaussEpp(valary_b.size());
//    valarray<double> valary_sol = gaussEpp.gauss_solver(A, valary_b);

// 调用求解方法来求解线性方程组
    SolveLinearEquations solveLinearEquations(A, valary_b, valary_b.size());
    DirectMethod directMethod(A, valary_b, valary_b.size());
    IterativeMethod iterativeMethod(A, valary_b, valary_b.size(), 1.0e-6, 1000);

    valarray<double> valary_sol;
    if (basis_type_try == 101 and basis_type_test == 101) {  // 1阶单元2节点用TMDA，1W单元用时3814.2ms，Gauss是3.27557e+06ms
        valary_sol = directMethod.tdma_solver();
    } else if (basis_type_try == 102 and basis_type_test == 102 and fast_method) {  // 1阶单元3节点改用迭代法，Gmres或BicgStab
        valary_sol = iterativeMethod.bicgstab_solver();  //81921.871000ms，1000高阶单元
    } else if (basis_type_try == 102 and basis_type_test == 102 and !fast_method) {
        valary_sol = iterativeMethod.gmres_solver();  //159257ms，1000高阶单元
    } else {
        cout << "103 not write" << endl;
        exit(1);
    }

    //  输出特征值/向量控制
    if (print_eigen) {
        int max_iter = 300;
        double eps = 1e-9;
        SolveEigenValue solveEigenValue(A, max_iter, eps);
//    myPrint.print_matrix(A);
        if (entire) {  // 完整输出 最大特征值和对应的特征向量，以及特征值矩阵和特征向量矩阵
            solveEigenValue.power_method();
            solveEigenValue.jacobi_rush();
        } else {  // 非完整输出 最大特征值，和对应的特征向量
            solveEigenValue.power_method();
        }
    }

//  输出计算结果控制，每个节点上都有结果
    if (fem_result) {
        cout << "***------Solution------***" << endl;
//    myPrint.print_matrix(vtr_sol);
        myPrint.print_vector(valary_sol);
    }

// 测试计算误差，调用精确解对比
    vector<vector<double>> vtr_sol = valary1d_to_vtr2d_trans1d(valary_sol);
    int n_basis = N_basis + 1;
    ErrorCompareAnalysis errorCompareAnalysis(n_basis, P);
    vector<double> error_list(0, 0);
    errorCompareAnalysis.get_entire_error_1d(vtr_sol, error_list, error_ary);

}

void FeSolverPoisson1D::fe_solver_poisson_2d(int &accuracy, bool fem_result = false) {

    MyPrint myPrint(accuracy);

    int N_mesh_v = int((right - left) / h[0]);
    int N_mesh_h = int((top - bottom) / h[1]);

    int N1_basis;
    int N2_basis;
    if (basis_type_try and basis_type_test == 201) {
        N1_basis = N_mesh_v;
        N2_basis = N_mesh_h;
    } else if (basis_type_try and basis_type_test == 202) {
        N1_basis = N_mesh_v * 2;
        N2_basis = N_mesh_h * 2;
    } else {
        cout << "Not write 203!" << endl;
    }


    vector<vector<double>> P;
    vector<vector<double>> T;

    GenPbTb genPbTb(left, right, bottom, top, h, basis_type_try, basis_type_test);
    genPbTb.gen_pb_tb_2d(P, T);

    size_t number_of_element = T[0].size();
    size_t N_basis = P[0].size() - 1;


    cout << "--------Global--------" << endl;
    cout << "Number of mesh = " << N_mesh_h * N_mesh_v * 2 << endl;
    cout << "Number of elements = " << number_of_element << endl;
    cout << "Number of nodes = " << N_basis+1 << endl;
    cout << "--------------------------------------------------------" << endl;


    int number_of_local_basis_fun_try = 0;
    vector<vector<double>> Tb_try;
    if (basis_type_try == 201) {
//    Pb_try = P
        Tb_try = T;
        number_of_local_basis_fun_try = 3;
    } else if (basis_type_try == 202) {
//    Pb_try = P
        Tb_try = T;
        number_of_local_basis_fun_try = 6;
    } else {
        cout << "pass" << endl;
    }

    int number_of_local_basis_fun_test = 0;
    vector<vector<double>> Tb_test;
    if (basis_type_test == 201) {
//    Pb_test = P
        Tb_test = T;
        number_of_local_basis_fun_test = 3;
    } else if (basis_type_test == 202) {
//    Pb_test = P
        Tb_test = T;
        number_of_local_basis_fun_test = 6;
    } else {
        cout << "pass" << endl;
    }

    vector<int> matrix_size = {int(N_basis) + 1, int(N_basis) + 1};
    int vector_size = int(N_basis) + 1;


    int basis_der_x_try = 1;
    int basis_der_y_try = 0;
    int basis_der_x_test = 1;
    int basis_der_y_test = 0;

    int basis_der_x_b_test = 0;
    int basis_der_y_b_test = 0;

    int number_of_gauss_point = int(number_of_integral);

    string str_test = "NiuLiuRui";

    AsmMtx1D asmMtxA1(str_test, matrix_size, int(number_of_element),
                      number_of_local_basis_fun_try, number_of_local_basis_fun_test,
                      P, T, Tb_test, Tb_try,
                      basis_type_try, basis_type_test,
                      basis_der_x_try, basis_der_x_test,
                      basis_der_y_try, basis_der_y_test,
                      number_of_gauss_point);
    vector<vector<double>> A1 = asmMtxA1.assembly_matrix_2d();

    basis_der_x_try = 0;
    basis_der_y_try = 1;
    basis_der_x_test = 0;
    basis_der_y_test = 1;
    AsmMtx1D asmMtxA2(str_test, matrix_size, int(number_of_element),
                      number_of_local_basis_fun_try, number_of_local_basis_fun_test,
                      P, T, Tb_test, Tb_try,
                      basis_type_try, basis_type_test,
                      basis_der_x_try, basis_der_x_test,
                      basis_der_y_try, basis_der_y_test,
                      number_of_gauss_point);
    vector<vector<double>> A2 = asmMtxA2.assembly_matrix_2d();
    vector<vector<double>> A;
    A = A1 + A2;

    AsmVtr1D asmVtr1D(str_test, vector_size, int(number_of_element), P, T, Tb_test,
                      number_of_local_basis_fun_test, basis_type_test, basis_der_x_b_test, basis_der_y_b_test,
                      number_of_gauss_point);

    vector<vector<double>> b = asmVtr1D.assembly_vector_2d();

    // 调用生成边界节点信息矩阵函数，传实参【单元个数整型】
    vector<vector<double>> boundary_nodes(0, vector<double>(0, 0));
    vector<vector<double>> boundary_edges(0, vector<double>(0, 0));
    BoundaryInformation boundaryInformation(N1_basis, N2_basis, N_mesh_v, N_mesh_h,
                                            boundary_nodes, boundary_edges);
    boundaryInformation.generate_boundary_nodes_2d();

// A, b, P放在函数方便引用改变形参
// 调用处理边界条件函数，传实参【未修改的整体刚度矩阵，未修改的载荷向量矩阵，边界节点信息矩阵】
    boundaryInformation.treat_boundary_condition(A, b, P, 2);

//    myPrint.print_matrix(A);
//    myPrint.print_matrix(b);

    valarray<double> valary_b = vtr2d_to_valary1d(b);
    SolveLinearEquations solveLinearEquations(A, valary_b, valary_b.size());
    IterativeMethod iterativeMethod(A, valary_b, valary_b.size(), 1.0e-6, 1000);

    valarray<double> valary_sol;
    valary_sol = iterativeMethod.bicgstab_solver();

//  输出计算结果控制，每个节点上都有结果
    cout << "***------Solution------***" << endl;
    if (fem_result) myPrint.print_vector(valary_sol);

    vector<vector<double>> vtr_sol = valary1d_to_vtr2d_trans1d(valary_sol);
    int n_basis = int(N_basis) + 1;
    ErrorCompareAnalysis errorCompareAnalysis(n_basis, P);
    errorCompareAnalysis.get_entire_error_2d(vtr_sol);

}

int main() {

//    struct timeval t1{}, t2{};
//    gettimeofday(&t1, nullptr);
//
//    int init_left = 0;
//    int init_right = 1;
//    vector<double> init_h = {0.25};
//    int init_try_type = 101;
//    int init_test_type = 101;
//    int init_integral_nun = 4;
//    int out_accuracy = 8;
//
//    FeSolverPoisson1D feSolverPoisson1D(init_left, init_right, init_h,
//                                        init_try_type, init_test_type,
//                                        init_integral_nun);
//
//    feSolverPoisson1D.fe_solver_poisson_1d(out_accuracy, false, false,
//                                           false, false, false, false);
//
//    gettimeofday(&t2, nullptr);
//    double time_cost = ((t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0);
//    std::cout << "time cost " << time_cost << "ms" << std::endl;


    struct timeval t1{}, t2{};
    gettimeofday(&t1, nullptr);

    int init_left = -1;
    int init_right = 1;
    int init_bottom = -1;
    int init_top = 1;
    vector<double> init_h = {1, 1};
    int init_try_type = 201;
    int init_test_type = 201;
    int init_integral_nun = 3;
    int out_accuracy = 8;

    FeSolverPoisson1D feSolverPoisson1D(init_left, init_right, init_bottom, init_top, init_h,
                                        init_try_type, init_test_type,
                                        init_integral_nun);

    feSolverPoisson1D.fe_solver_poisson_2d(out_accuracy, false);

    gettimeofday(&t2, nullptr);
    double time_cost = ((t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0);
    std::cout << "time cost " << time_cost << "ms" << std::endl;


    return 0;
}
//10000 Elements
// Gauss vs Tdma
//1.03431e-08%
//time cost 3.27557e+06ms

//1.74058e-09%
//time cost 3814.2ms