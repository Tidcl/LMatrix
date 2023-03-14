//
// Created by admin on 2023/1/5.
//

#ifndef MAIN_CPP_FEBASFUN1D_H
#define MAIN_CPP_FEBASFUN1D_H

#include <iostream>
#include <vector>

using namespace std;

class BasisFun1D {
private:
    double x;
    double y;
    int basis_type;
    int basis_index;
    int basis_der_x;
    int basis_der_y;
    vector<vector<double>> vertices;

public:
    BasisFun1D(double &x, vector<vector<double>> &vertices, int &basis_type, int &basis_index, int &basis_der_x);

    BasisFun1D(double &x, double &y,vector<vector<double>> &vertices,
               int &basis_type, int &basis_index, int &basis_der_x, int &basis_der_y);

    ~BasisFun1D() = default;

    double fe_basis_fun_1d();

    double fe_reference_basis_fun_2d(double xh, double yh,int ref_basis_type, int ref_basis_index,
                                     int ref_basis_der_x, int ref_basis_der_y) const;

    double fe_local_basis_fun_2d();
};


#endif //MAIN_CPP_FEBASFUN1D_H
