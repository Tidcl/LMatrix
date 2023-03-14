//
// Created by admin on 2023/1/6.
//

#include <cmath>
#include "../Tools/PrintTools.h"
#include "BoundaryInformation.h"

// 定义了边界节点信息矩阵函数，传参【单元个数 == 整体基函数个数】
BoundaryInformation::BoundaryInformation(int &N_basis, vector<vector<double>> &boundary_nodes) {
    this->N_basis = N_basis;
    this->boundary_nodes = boundary_nodes;
//    this->A = A;
//    this->b = b;
//    this->P = P;
}

BoundaryInformation::BoundaryInformation(int &N1_basis, int &N2_basis, int &N1_partition, int &N2_partition,
                                         vector<vector<double>> &boundary_nodes,
                                         vector<vector<double>> &boundary_edges) {
    this->N1_basis = N1_basis;
    this->N2_basis = N2_basis;
    this->N1_partition = N1_partition;
    this->N2_partition = N2_partition;
    this->boundary_nodes = boundary_nodes;
    this->boundary_edges = boundary_edges;
}

void BoundaryInformation::generate_boundary_nodes() {
//# 生成boundary_nodes矩阵的数据结构
    boundary_nodes = vector<vector<double>>(3, vector<double>(2, 0));
// 类型(狄利克雷或纽曼)
    boundary_nodes[0][0] = -1;
    boundary_nodes[0][1] = -1;
// 边界节点的全局索引
    boundary_nodes[1][0] = 0;
    boundary_nodes[1][1] = N_basis;
// 边界节点的法线方向
    boundary_nodes[2][0] = -1;
    boundary_nodes[2][1] = 1;

}

void BoundaryInformation::generate_boundary_nodes_2d() {
    int nbn = 2 * (N1_basis + N2_basis);
//    boundary_nodes = np.zeros([2, nbn])
    boundary_nodes = vector<vector<double>>(2, vector<double>(nbn, 0));
    for (int i = 0; i < boundary_nodes[0].size(); ++i) {
        boundary_nodes[0][i] = -1;
    }

    for (int k = 1; k < N1_basis + 1; ++k) {
        boundary_nodes[1][k - 1] = (k - 1) * (N2_basis + 1) + 1 - 1;
    }

    for (int k = N1_basis + 1; k < N1_basis + N2_basis + 1; ++k) {
        boundary_nodes[1][k - 1] = N1_basis * (N2_basis + 1) + k - N1_basis - 1;
    }

    for (int k = N1_basis + N2_basis + 1; k < 2 * N1_basis + N2_basis + 1; ++k) {
        boundary_nodes[1][k - 1] = (2 * N1_basis + N2_basis + 2 - k) * (N2_basis + 1) - 1;
    }

    for (int k = 2 * N1_basis + N2_basis + 1; k < nbn + 1; ++k) {
        boundary_nodes[1][k - 1] = 2 * N1_basis + 2 * N2_basis + 2 - k - 1;
    }
}

void BoundaryInformation::generate_boundary_edges_2d() {
    int nbe = 2 * (N1_partition + N2_partition);
//    boundary_edges = np.zeros([4, nbe])
    boundary_edges = vector<vector<double>>(4, vector<double>(nbe, 0));
    for (int i = 0; i < boundary_edges[0].size(); ++i) {
        boundary_edges[0][i] = -1;
    }

    for (int k = 1; k < N1_partition + 1; ++k) {
        boundary_edges[1][k - 1] = (k - 1) * 2 * N2_partition + 1 - 1;
        boundary_edges[2][k - 1] = (k - 1) * (N2_partition + 1) + 1 - 1;
        boundary_edges[3][k - 1] = k * (N2_partition + 1) + 1 - 1;
    }

    for (int k = N1_partition + 1; k < N1_partition + N2_partition + 1; ++k) {
        boundary_edges[1][k - 1] = (N1_partition - 1) * 2 * N2_partition + 2 * (k - N1_partition) - 1;
        boundary_edges[2][k - 1] = N1_partition * (N2_partition + 1) + k - N1_partition - 1;
        boundary_edges[3][k - 1] = N1_partition * (N2_partition + 1) + k - N1_partition + 1 - 1;
    }

    for (int k = N1_partition + N2_partition + 1; k < 2 * N1_partition + N2_partition + 1; ++k) {
        boundary_edges[1][k - 1] = (2 * N1_partition + N2_partition + 1 - k) * 2 * N2_partition - 1;
        boundary_edges[2][k - 1] = (2 * N1_partition + N2_partition + 2 - k) * (N2_partition + 1) - 1;
        boundary_edges[3][k - 1] = (2 * N1_partition + N2_partition + 1 - k) * (N2_partition + 1) - 1;
    }

    for (int k = 2 * N1_partition + N2_partition + 1; k < nbe + 1; ++k) {
        boundary_edges[1][k - 1] = 2 * (2 * N1_partition + 2 * N2_partition + 1 - k) - 1 - 1;
        boundary_edges[2][k - 1] = 2 * N1_partition + 2 * N2_partition + 2 - k - 1;
        boundary_edges[3][k - 1] = 2 * N1_partition + 2 * N2_partition + 1 - k - 1;
    }
}

// 处理边界，传参【A矩阵，b矩阵，boundary_nodes矩阵】
void BoundaryInformation::treat_boundary_condition(vector<vector<double>> &A, vector<vector<double>> &b,
                                                   vector<vector<double>> &P, int dim_mark) {
    char temp = 'b';
    double init_fx = 0;
//    获取boundary_nodes矩阵的列数，nbn = 2
    int nbn = int(boundary_nodes[0].size());

    for (int k = 0; k < nbn; ++k) {
        if (boundary_nodes[0][k] == -1) {  // 处理狄利克雷边界
            int i = int(boundary_nodes[1][k]);
//            A[i][:] = 0;
            for (int j = 0; j < A[0].size(); ++j) {  // A对应行变0
                A[i][j] = 0;
            }
            A[i][i] = 1;  // A行列对应值变1
            if (dim_mark == 1){
                b[i][0] = bc_function(temp, int(P[0][i]), init_fx);  // # b对应行变成边界给的值，调用边界函数，传实参【】
            } else if (dim_mark ==2){
                b[i][0] = bc_function(temp, double(P[0][i]),  double(P[1][i]),init_fx);
            }else{
                cout<<"pass"<<endl;
            }
        }
    }
//    int accuracy = 6;
//    MyPrint myPrint(accuracy);
//    myPrint.print_matrix(A);
//    cout<<&A<<"nei A "<<endl;
//    myPrint.print_matrix(b);
//    myPrint.print_matrix(P);
//    print_matrix_bc(A);
//    print_matrix_bc(b);
}

// 定义了边界函数，定义形参【边界选择标记】上面的调用基函数 用类参数和方法参数传参，这里调用同时传参即可
double BoundaryInformation::bc_function(char &bc_select_mark, int x, double &fx) {
    if (bc_select_mark == 'b') {
        int left = 0;
        int right = 1;
        if (x == left) {
            fx = 0;  // 人为给的
        } else if (x == right) {
            fx = cos(1);  // 人为给的
        } else {
            fx = 110;
            cout << "That's what I did on purpose! 110." << endl;
        }
    } else {
        fx = 321;
        cout << "That's what I did on purpose! 321." << endl;
    }
    return fx;
}

double BoundaryInformation::bc_function(char &bc_select_mark, double x, double y, double &fx) {
    if (bc_select_mark == 'b') {
        if (x == -1) {
            fx = -1.5 * y * (1 - y) * exp(-1 + y);  // 人为给的
        } else if (x == 1) {
            fx = 0.5 * y * (1 - y) * exp(1 + y);  // 人为给的
        } else if (y == -1) {
            fx = -2 * x * (1 - x / 2.0) * exp(x - 1);
        } else if (y == 1) {
            fx = 0;
        }
    } else {
        fx = 321;
        cout << "That's what I did on purpose! 321." << endl;
    }
    return fx;
}

//int main() {
//    //test 4 pass
//    int N_basis = 4;
//    vector<vector<double>> boundary_nodes(0, vector<double>(0, 0));
//    vector<vector<double>> A = {{ 4.54440667, -4.54440667,  0.0       ,  0.0       ,  0.0       },
//                                {-4.54440667, 10.37954033, -5.83513366,  0.0       ,  0.0       },
//                                { 0.0       , -5.83513366, 13.3275936 , -7.49245993,  0.0       },
//                                { 0.0       ,  0.0       , -7.49245993, 17.11296892, -9.62050899},
//                                { 0.0       ,  0.0       ,  0.0       , -9.62050899,  9.62050899}};
//    vector<vector<double>> b{{-0.09856355},
//                             {-0.04020103},
//                             { 0.33110004},
//                             { 0.91526931},
//                             { 0.71105657}};
//    vector<vector<double>> P{{0.0, 0.25, 0.5, 0.75, 1.0}};
//
//    BoundaryInformation boundaryInformation(N_basis, boundary_nodes, A, b, P);
//    boundaryInformation.generate_boundary_nodes();
//    boundaryInformation.treat_boundary_condition();
//}