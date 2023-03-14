//
// Created by admin on 2023/1/5.
//

#include <cmath>
#include "FeBasisFunction.h"

using namespace std;

// 定义了生成局部基函数的函数，定义形参【单元上的任意坐标x，单元节点坐标，基函数的在节点上的类型，基函数的在节点上的索引，基函数的阶数】
BasisFun1D::BasisFun1D(double &x, vector<vector<double>> &vertices,
                       int &basis_type, int &basis_index, int &basis_der_x) {
    this->x = x;
    this->vertices = vertices;
    this->basis_type = basis_type;
    this->basis_index = basis_index;
    this->basis_der_x = basis_der_x;
}

double BasisFun1D::fe_basis_fun_1d() {
    double result = 0.0;
//    double h = vertices[0][-1] - vertices[0][0];
/*
    h = vertices[0][1] - vertices[0][0]  ==> h = 0.25
    这里是天坑中的天坑，基函数是有限元概念，即1个单元，2个子单元，3个节点，3个局部基函数
    所有当1阶单元时，vertices每一行存了每一个单元的2个坐标，但是2阶单元时，存了每一行存了3个坐标
    同时顺序是生成T矩阵时的左中右，因此只有最左和最右是其真正的单元节点坐标
    h是代表真正意义的每一个单元的尺寸！
 */
    double h = vertices[0].back() - vertices[0].front();

    if (basis_type == 101) {
        if (basis_index == 0) {
            if (basis_der_x == 0) {
                result = (vertices[0][1] - x) / h;
            } else if (basis_der_x == 1) {
                result = -1 / h;
            } else if (basis_der_x >= 2) {
                result = 0;
            } else {
                cout << "Error der order" << endl;
            }
        } else if (basis_index == 1)
            if (basis_der_x == 0) {
                result = (x - vertices[0][0]) / h;
            } else if (basis_der_x == 1) {
                result = 1 / h;
            } else if (basis_der_x >= 2) {
                result = 0;
            } else {
                cout << "Error der order" << endl;
            }
        else {
            cout << "Error!" << endl;
        }
    } else if (basis_type == 102) {

        if (basis_index == 0) {
            if (basis_der_x == 0) {
                result = 2 * pow((x - vertices[0][0]) / h, 2) - 3 * ((x - vertices[0][0]) / h) + 1;
            } else if (basis_der_x == 1) {
                result = 4 / pow(h, 2) * (x - vertices[0][0]) - 3 / h;
            } else if (basis_der_x == 2) {
                result = 4 / pow(h, 2);
            } else if (basis_der_x > 2) {
                result = 0;
            } else {

            }
        } else if (basis_index == 1) {
            if (basis_der_x == 0) {
                result = -4 * pow((x - vertices[0][0]) / h, 2) + 4 * ((x - vertices[0][0]) / h);
            } else if (basis_der_x == 1) {
                result = -8 / pow(h, 2) * (x - vertices[0][0]) + 4 / h;
            } else if (basis_der_x == 2) {
                result = -8 / pow(h, 2);
            } else if (basis_der_x > 2) {
                result = 0;
            } else {
                cout << "Error der order" << endl;
            }
        } else if (basis_index == 2) {
            if (basis_der_x == 0) {
                result = 2 * pow((x - vertices[0][0]) / h, 2) - ((x - vertices[0][0]) / h);
            } else if (basis_der_x == 1) {
                result = 4 / pow(h, 2) * (x - vertices[0][0]) - 1 / h;
            } else if (basis_der_x == 2) {
                result = 4 / pow(h, 2);
            } else if (basis_der_x > 2) {
                result = 0;
            } else {
                cout << "Error der order" << endl;
            }
        }
    }
    return result;
}

BasisFun1D::BasisFun1D(double &x, double &y, vector<vector<double>> &vertices, int &basis_type, int &basis_index,
                       int &basis_der_x, int &basis_der_y) {
    this->x = x;
    this->y = y;
    this->vertices = vertices;
    this->basis_type = basis_type;
    this->basis_index = basis_index;
    this->basis_der_x = basis_der_x;
    this->basis_der_y = basis_der_y;
}

double BasisFun1D::fe_reference_basis_fun_2d(double xh, double yh, int ref_basis_type, int ref_basis_index,
                                             int ref_basis_der_x, int ref_basis_der_y) const {

    double result = 0.0;

    if (ref_basis_type == 201) {
        if (ref_basis_index == 0) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = 1 - xh - yh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = -1;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = -1;
        } else if (ref_basis_index == 1) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = xh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = 1;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = 0;

        } else if (ref_basis_index == 2) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = yh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = 0;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = 1;
        } else {
            cout << "not write" << endl;
        }
    } else if (ref_basis_type == 202) {
        if (ref_basis_index == 0) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = 2 * pow(xh, 2) + 2 * pow(yh, 2) + 4 * xh * yh - 3 * xh - 3 * yh + 1;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = 4 * xh + 4 * yh - 3;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = 4 * yh + 4 * xh - 3;
            if (ref_basis_der_x == 2 and ref_basis_der_y == 0)
                result = 4;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 2)
                result = 4;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 1)
                result = 4;
        } else if (ref_basis_index == 1) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = 2 * pow(xh, 2) - xh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = 4 * xh - 1;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = 0;
            if (ref_basis_der_x == 2 and ref_basis_der_y == 0)
                result = 4;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 2)
                result = 0;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 1)
                result = 0;
        } else if (basis_index == 2) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = 2 * pow(yh, 2) - yh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = 0;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = 4 * yh - 1;
            if (ref_basis_der_x == 2 and ref_basis_der_y == 0)
                result = 0;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 2)
                result = 4;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 1)
                result = 0;
        } else if (basis_index == 3) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = 4 * xh - 4 * pow(xh, 2) - 4 * xh * yh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = 4 - 8 * xh - 4 * yh;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = -4 * xh;
            if (ref_basis_der_x == 2 and ref_basis_der_y == 0)
                result = -8;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 2)
                result = 0;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 1)
                result = -4;
        } else if (basis_index == 4) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = 4 * xh * yh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = 4 * yh;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = 4 * xh;
            if (ref_basis_der_x == 2 and ref_basis_der_y == 0)
                result = 0;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 2)
                result = 0;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 1)
                result = 4;
        } else if (basis_index == 5) {
            if (ref_basis_der_x == 0 and ref_basis_der_y == 0)
                result = 4 * yh - 4 * pow(yh, 2) - 4 * xh * yh;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 0)
                result = -4 * yh;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 1)
                result = 4 - 8 * yh - 4 * xh;
            if (ref_basis_der_x == 2 and ref_basis_der_y == 0)
                result = 0;
            if (ref_basis_der_x == 0 and ref_basis_der_y == 2)
                result = -8;
            if (ref_basis_der_x == 1 and ref_basis_der_y == 1)
                result = -4;
        } else {
            cout << "not write" << endl;
        }
    } else {
        cout << "not write" << endl;
    }
    return result;
}

double BasisFun1D::fe_local_basis_fun_2d() {
    double J_11 = vertices[0][1] - vertices[0][0];
    double J_12 = vertices[0][2] - vertices[0][0];
    double J_21 = vertices[1][1] - vertices[1][0];
    double J_22 = vertices[1][2] - vertices[1][0];

    double jacobi_det = J_11 * J_22 - J_12 * J_21;
    double xh = (J_22 * (x - vertices[0][0]) - J_12 * (y - vertices[1][0])) / jacobi_det;
    double yh = (-J_21 * (x - vertices[0][0]) + J_11 * (y - vertices[1][0])) / jacobi_det;

    double result = 0.0;

    if (basis_der_x == 0 and basis_der_y == 0)
        result = fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, basis_der_x, basis_der_y);
    else if (basis_der_x == 1 and basis_der_y == 0)
        result = (J_22 / jacobi_det) * fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 1, 0) + \
                   (-J_21 / jacobi_det) * fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 0, 1);
    else if (basis_der_x == 0 and basis_der_y == 1)
        result = (-J_12 / jacobi_det) * fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 1, 0) + \
                   (J_11 / jacobi_det) * fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 0, 1);
    else if (basis_der_x == 2 and basis_der_y == 0)
        result =
                (pow(J_22, 2) / pow(jacobi_det, 2)) * fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 2, 0) + \
                   (pow(J_21, 2) / pow(jacobi_det, 2)) *
                   fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 0, 2) + \
                   (-2 * J_21 * J_22 / pow(jacobi_det, 2)) *
                   fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 1, 1);
    else if (basis_der_x == 0 and basis_der_y == 2)
        result =
                (pow(J_12, 2) / pow(jacobi_det, 2)) * fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 2, 0) + \
                   (pow(J_11, 2) / pow(jacobi_det, 2)) *
                   fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 0, 2) + \
                   (-2 * J_11 * J_12 / pow(jacobi_det, 2)) *
                   fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 1, 1);
    else if (basis_der_x == 1 and basis_der_y == 1)
        result =
                (-J_22 * J_12 / pow(jacobi_det, 2)) * fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 2, 0) + \
                   (-J_21 * J_11 / pow(jacobi_det, 2)) *
                   fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 0, 2) + \
                   (J_21 * J_12 + J_11 * J_22 / pow(jacobi_det, 2)) *
                   fe_reference_basis_fun_2d(xh, yh, basis_type, basis_index, 1, 1);
    else {
        cout << "not write" << endl;
    }
    return result;
}