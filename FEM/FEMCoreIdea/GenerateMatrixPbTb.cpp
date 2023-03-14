//
// Created by Administrater on 2023/1/4.
//
#include <cmath>
#include <vector>
#include "GenerateMatrixPbTb.h"
#include "../Tools/PrintTools.h"

// 定义生成网格信息矩阵函数，定义形参【左侧点坐标（1个），右侧点坐标（1个），尺寸，试探基(形)函数类型，测试基(形)函数类型】
GenPbTb::GenPbTb(double left, double right, double h, int basis_type_try, int basis_type_test) {
    this->left = left;
    this->right = right;
    this->h = h;
    this->basis_type_try = basis_type_try;
    this->basis_type_test = basis_type_test;
}

GenPbTb::GenPbTb(double left, double right, double bottom, double top, vector<double> &h2, int basis_type_try,
                 int basis_type_test) {
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->h2 = h2;
    this->basis_type_try = basis_type_try;
    this->basis_type_test = basis_type_test;
}

void GenPbTb::gen_pb_tb_1d(vector<vector<double>> &P, vector<vector<double>> &T, double &h_basis) const {
//  通过四舍五入判断余数，确定网格数/单元数  N = 4
    int N;
    double remainder;
    int big_number = 1000000;
    remainder = fmod((right - left) * big_number, h * big_number);
    if (remainder == 0) {
        N = int(((right - left) * big_number) / (h * big_number));
    } else {
        N = int(((right - left) * big_number) / (h * big_number));
        if (remainder >= 0.5 * big_number) {
            N = N + 1;
        } else {
            cout << "Error: I haven't written it yet!" << endl;
        }
    }

//  修改P，T矩阵数据，这里开始，若是不能整除，执行了四舍五入，应该修改尺寸
    if (basis_type_try and basis_type_test == 101) {
//        vector<double> P(N + 1, 0);
//        vector<vector<double>> T(2, vector<double>(N, 0));
        P = vector<vector<double>>(1, vector<double>(N + 1, 0));
        T = vector<vector<double>>(2, vector<double>(N, 0));
//        print_matrix(P);
//        print_matrix(T);
        for (int i = 0; i < N + 1; ++i) {
            P[0][i] = left + i * h;
        }
        for (int j = 0; j < N; ++j) {
            T[0][j] = j + 0;
            T[1][j] = j + 1;
        }
        h_basis = h;
//        print_matrix(P);
//        print_matrix(T);
//        cout<<h_basis<<endl;
    } else if (basis_type_try and basis_type_test == 102) {
        P = vector<vector<double>>(1, vector<double>(2 * N + 1, 0));
        T = vector<vector<double>>(3, vector<double>(N, 0));
//        print_matrix(P);
//        print_matrix(T);
        for (int i = 0; i < 2 * N + 1; ++i) {
            P[0][i] = left + i * h / 2.0;
        }
        for (int j = 0; j < N; ++j) {
            T[0][j] = 2 * j;
            T[1][j] = 2 * j + 1;
            T[2][j] = 2 * (j + 1);
        }
        h_basis = h / 2;
//        print_vector(P[0]);
//        print_matrix(T);
//        cout<<h_basis<<endl;
    } else {
        cout << "Error: I haven't written it yet!" << endl;
        exit(1);
    }
}

void GenPbTb::gen_pb_tb_2d(vector<vector<double>> &P, vector<vector<double>> &T) {
    int col, row;
    vector<vector<double>> Q;
    double N1 = (right - left) / h2[0];
    double N2 = (top - bottom) / h2[1];
    if (basis_type_try == 201 and basis_type_test == 201) {
        double node_num = (N1 + 1) * (N2 + 1);

        P = vector<vector<double>>(2, vector<double>(int(node_num), 0));
        T = vector<vector<double>>(3, vector<double>(int(2 * N1 * N2), 0));
        Q = vector<vector<double>>(int(N1 + 1), vector<double>(int(N2 + 1), 0));

        for (int j = 0; j < node_num; ++j) {
            if ((j + 1) % int(N2 + 1) == 0) {
                P[0][j] = left + ((j + 1) / (N2 + 1) - 1) * h2[0];
                P[1][j] = top;
            } else {
                P[0][j] = left + ((j + 1) / int(N2 + 1)) * h2[0];
                P[1][j] = bottom + ((j + 1) % int(N2 + 1) - 1) * h2[1];
            }
        }
        for (int i = 0; i < int(N1 + 1); ++i) {
            for (int j = 0; j < int(N2 + 1); ++j) {
                Q[i][j] = int(i * (N2 + 1) + j + 1);
            }
        }

        for (int n = 0; n < int(N1 * N2); ++n) {
            if ((n + 1) % int(N2) == 0) {
                row = int(N2) - 1;
                col = int((n + 1) / int(N2)) - 1;
            } else {
                row = int((n + 1) % int(N2)) - 1;
                col = int((n + 1) / int(N2)) + 1 - 1;
            }

            T[0][2 * n] = Q[col][row] - 1;
            T[1][2 * n] = Q[col + 1][row] - 1;
            T[2][2 * n] = Q[col][row + 1] - 1;

            T[0][2 * n + 1] = Q[col][row + 1] - 1;
            T[1][2 * n + 1] = Q[col + 1][row] - 1;
            T[2][2 * n + 1] = Q[col + 1][row + 1] - 1;
        }
//        int sds = 2;
//        MyPrint myPrint(sds);
//        myPrint.print_matrix(P);
//        myPrint.print_matrix(T);
//        myPrint.print_matrix(Q);
    } else if (basis_type_try == 202 and basis_type_test == 202) {

        vector<double> dh = {h2[0] / 2, h2[1] / 2};
        double dN1 = N1 * 2;
        double dN2 = N2 * 2;

        double node_num = (dN1 + 1) * (dN2 + 1);

        P = vector<vector<double>>(2, vector<double>(int(node_num), 0));
        T = vector<vector<double>>(6, vector<double>(int(2 * N1 * N2), 0));
        Q = vector<vector<double>>(int(dN1 + 1), vector<double>(int(dN2 + 1), 0));

        for (int j = 0; j < node_num; ++j) {
            if ((j + 1) % int(dN2 + 1) == 0) {
                P[0][j] = left + ((j + 1) / (dN2 + 1) - 1) * dh[0];
                P[1][j] = top;
            } else {
                P[0][j] = left + (j + 1) / int(dN2 + 1) * dh[0];
                P[1][j] = bottom + ((j + 1) % int(dN2 + 1) - 1) * dh[1];
            }
        }

        for (int i = 0; i < int(dN1 + 1); ++i) {
            for (int j = 0; j < int(dN2 + 1); ++j) {
                Q[i][j] = i * (dN2 + 1) + j + 1;
            }
        }

        for (int n = 0; n < int(N1 * N2); ++n) {
            if ((n + 1) % int(N2) == 0) {
                row = int(N2) - 1;
                col = int((n + 1) / N2) - 1;
            } else {
                row = int((n + 1) % int(N2)) - 1;
                col = int((n + 1) / N2) + 1 - 1;
            }
            T[0][2 * n] = Q[2 * col - 1 + 1][2 * row - 1 + 1] - 1;
            T[1][2 * n] = Q[2 * col + 1 + 1][2 * row - 1 + 1] - 1;
            T[2][2 * n] = Q[2 * col - 1 + 1][2 * row + 1 + 1] - 1;
            T[3][2 * n] = Q[2 * col + 0 + 1][2 * row - 1 + 1] - 1;
            T[4][2 * n] = Q[2 * col + 0 + 1][2 * row + 0 + 1] - 1;
            T[5][2 * n] = Q[2 * col - 1 + 1][2 * row + 0 + 1] - 1;

            T[0][2 * n + 1] = Q[2 * col - 1 + 1][2 * row + 1 + 1] - 1;
            T[1][2 * n + 1] = Q[2 * col + 1 + 1][2 * row - 1 + 1] - 1;
            T[2][2 * n + 1] = Q[2 * col + 1 + 1][2 * row + 1 + 1] - 1;
            T[3][2 * n + 1] = Q[2 * col + 0 + 1][2 * row + 0 + 1] - 1;
            T[4][2 * n + 1] = Q[2 * col + 1 + 1][2 * row + 0 + 1] - 1;
            T[5][2 * n + 1] = Q[2 * col + 0 + 1][2 * row + 1 + 1] - 1;
        }
//        int sds = 2;
//        MyPrint myPrint(sds);
//        myPrint.print_matrix(P);
//        myPrint.print_matrix(T);
//        myPrint.print_matrix(Q);
    } else {
        cout << "not write!" << endl;
    }
}

//int main() {
//    //test 1 pass
//    double val_left = -1;
//    double val_right = 1;
//    double val_bottom = -1;
//    double val_top = 1;
//    vector<double> val_h = {1, 1};
//    int vak_type_try = 201;
//    int vak_type_test = 201;
//
//    GenPbTb test_fun_1(val_left, val_right, val_bottom, val_top,
//                       val_h, vak_type_try, vak_type_test);
//
//    vector<vector<double>> P(0, vector<double>(0, 0));
//    vector<vector<double>> T(0, vector<double>(0, 0));
//    double h_basis = 0;
//    test_fun_1.gen_pb_tb_2d(P, T);
//
//    cout << "hahaha" << endl;
//    return 0;
//}