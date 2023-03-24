//
// Created by admin on 2023/1/6.
//

#ifndef MAIN_CPP_BDYINFOR_H
#define MAIN_CPP_BDYINFOR_H

#include <iostream>
#include "../NMatrixDefine.h"

using namespace std;


class BoundaryInformation {
private:
    int N_basis;
    matrixv boundary_nodes;

    int N1_basis;
    int N2_basis;
    int N1_partition;
    int N2_partition;
    matrixv boundary_edges;
//    函数引用传参不改变地址，类引用传参因为构造的原因，this->A = A，会改变地址
//    matrixv A;
//    matrixv b;
//    matrixv P;

public:
//    BoundaryInformation(int &N_basis, matrixv &boundary_nodes,
//                        matrixv &A, matrixv &b,
//                        matrixv &P);

    BoundaryInformation(int &N_basis, matrixv &boundary_nodes);

    BoundaryInformation(int &N1_basis, int &N2_basis, int &N1_partition, int &N2_partition,
                        matrixv &boundary_nodes, matrixv &boundary_edges);

    ~BoundaryInformation() = default;

    void generate_boundary_nodes();

    void generate_boundary_nodes_2d();

    void generate_boundary_edges_2d();

    void treat_boundary_condition(matrixv &A, matrixv &b,
                                  matrixv &P, int dim_mark);

    static double bc_function(char &bc_select_mark, int x, double &fx);

    static double bc_function(char &bc_select_mark, double x, double y,double &fx);
};

//void print_matrix_bc(matrixv &matrix) {
//    for (int i = 0; i < matrix.size(); ++i) {
//        for (int j = 0; j < matrix[0].size(); ++j) {
//            cout << matrix[i][j] << '\t' << '\t';
//        }
//        cout << endl;
//    }
//}


#endif //MAIN_CPP_BDYINFOR_H
