//
// Created by admin on 2023/1/6.
//

#ifndef MAIN_CPP_BDYINFOR_H
#define MAIN_CPP_BDYINFOR_H

#include <iostream>
#include <vector>

using namespace std;


class BoundaryInformation {
private:
    int N_basis;
    vector<vector<double>> boundary_nodes;

    int N1_basis;
    int N2_basis;
    int N1_partition;
    int N2_partition;
    vector<vector<double>>boundary_edges;
//    函数引用传参不改变地址，类引用传参因为构造的原因，this->A = A，会改变地址
//    vector<vector<double>> A;
//    vector<vector<double>> b;
//    vector<vector<double>> P;

public:
//    BoundaryInformation(int &N_basis, vector<vector<double>> &boundary_nodes,
//                        vector<vector<double>> &A, vector<vector<double>> &b,
//                        vector<vector<double>> &P);

    BoundaryInformation(int &N_basis, vector<vector<double>> &boundary_nodes);

    BoundaryInformation(int &N1_basis, int &N2_basis, int &N1_partition, int &N2_partition,
                        vector<vector<double>> &boundary_nodes, vector<vector<double>> &boundary_edges);

    ~BoundaryInformation() = default;

    void generate_boundary_nodes();

    void generate_boundary_nodes_2d();

    void generate_boundary_edges_2d();

    void treat_boundary_condition(vector<vector<double>> &A, vector<vector<double>> &b,
                                  vector<vector<double>> &P, int dim_mark);

    static double bc_function(char &bc_select_mark, int x, double &fx);

    static double bc_function(char &bc_select_mark, double x, double y,double &fx);
};

//void print_matrix_bc(vector<vector<double>> &matrix) {
//    for (int i = 0; i < matrix.size(); ++i) {
//        for (int j = 0; j < matrix[0].size(); ++j) {
//            cout << matrix[i][j] << '\t' << '\t';
//        }
//        cout << endl;
//    }
//}


#endif //MAIN_CPP_BDYINFOR_H
