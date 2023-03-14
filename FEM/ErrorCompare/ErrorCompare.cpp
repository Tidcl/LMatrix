//
// Created by admin on 2023/1/7.
//

#include <cmath>
#include "ErrorCompare.h"
#include "../Tools/PrintTools.h"


// 定义了误差计算函数，定义形参【数值解，整体基函数个数=单元个数，左侧端点，整个基函数尺寸=单元尺寸】
ErrorCompareAnalysis::ErrorCompareAnalysis(int &num_node, vector<vector<double>> &P) {
//    this->solution = solution;
    this->P = P;
    this->num_node = num_node;
//    this->h_basis = h_basis;
}


double ErrorCompareAnalysis::get_entire_error_1d(vector<vector<double>> &solution,
                                                 vector<double> &error_list,
                                                 bool show_error_list = false) const {
    int accuracy = 8;
    MyPrint myPrint(accuracy);
    double max_error = 0;
    error_list = vector<double>(0, 0.0);
    vector<double> exact_list(0, 0.0);
    vector<vector<double>> temp(0, vector<double>(solution[0].size(), 0.0));

    for (int i = 0; i < num_node; ++i) {
        double temp_error = solution[i][0] - exact_solution(P[0][i]);
        error_list.push_back(temp_error);
        exact_list.push_back(exact_solution(P[0][i]));
        if (fabs(max_error) < fabs(temp_error)) {
            max_error = temp_error * 100;
        }
    }
    if (show_error_list) {
        cout << "***----The MaxError----***" << endl;
        cout << fabs(max_error) << "%" << endl;
        cout << "***------Exact_list------***" << endl;
        myPrint.print_vector(exact_list);
        cout << "***------Error_list------***" << endl;
        myPrint.print_vector(error_list);
//        printVector(error_list);
//        printVector(exact_list);
    } else {
        cout << "***----The MaxError----***" << endl;
        cout << fabs(max_error) << "%" << endl;
    }
    return fabs(max_error);
}

double ErrorCompareAnalysis::get_entire_error_2d(vector<vector<double>> &solution) {
    double max_error = 0;
    for (int i = 0; i < num_node; ++i) {
        double temp_error = solution[i][0] - exact_solution(P[0][i], P[1][i]);
        if (fabs(max_error) < fabs(temp_error)) {
            max_error = temp_error * 100;
        }
    }
    cout << "***----The MaxError----***" << endl;
    cout << fabs(max_error) << "%" << endl;
    return fabs(max_error);
}

// 定义了精确解函数对比，定义形参【点坐标】
double ErrorCompareAnalysis::exact_solution(double x) {
    double exact_result = x * cos(x);
    return exact_result;
}

double ErrorCompareAnalysis::exact_solution(double x, double y) {
    double exact_result =x * y * (1 - x / 2) * (1 - y) * exp(x + y);
    return exact_result;
}

//int main() {
//    //test 6 pass
//    vector<vector<double>> solution = {{0.000000},
//                                       {0.244117},
//                                       {0.441125},
//                                       {0.550364},
//                                       {0.540302}};
//    int n_basis = 5;
//    int left = 0;
//    double h_basis = 0.25;
//
//    vector<double>error_list(0,0);
//
//    ErrorCompareAnalysis errorCompareAnalysis(solution, n_basis, left, h_basis);
//    errorCompareAnalysis.get_entire_error_1d(error_list, true);
//}
