//
// Created by Administrater on 2023/1/7.
//

#ifndef RUNSCRIPY_PY_FEMSOLVER1D_H
#define RUNSCRIPY_PY_FEMSOLVER1D_H

#include <iostream>
#include "../NMatrixDefine.h"

using namespace std;

class FeSolverPoisson1D {
public:
    FeSolverPoisson1D(int left, int right,int bottom, int top,vecd&h,
                      int basis_type_try, int basis_type_test,
                      int number_of_integral);

    FeSolverPoisson1D(int left, int right, vecd&h,
                      int basis_type_try, int basis_type_test,
                      int number_of_integral);

    ~FeSolverPoisson1D() = default;


    void fe_solver_poisson_1d(int &accuracy, bool print_A, bool print_eigen,
                              bool entire, bool fem_result, bool error_ary, bool fast_method);

    void fe_solver_poisson_2d(int &accuracy,bool fem_result);
private:
    double left, right, bottom{},top{};
    int number_of_integral;
    int basis_type_try, basis_type_test;
    vecd h;
};


#endif //RUNSCRIPY_PY_FEMSOLVER1D_H
