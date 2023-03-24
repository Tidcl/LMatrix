//
// Created by admin on 2023/1/7.
//

#ifndef MAIN_CPP_ERRORCOMPARE_H
#define MAIN_CPP_ERRORCOMPARE_H

#include <iostream>
#include "../NMatrixDefine.h"

using namespace std;

class ErrorCompareAnalysis {
private:
    matrixv P;
    int num_node;

public:
    ErrorCompareAnalysis(int &num_node, matrixv &P);

    ~ErrorCompareAnalysis() = default;

public:

    double get_entire_error_1d(matrixv &solution,
                               vecd &error_list, bool show_error_list) const;

    double get_entire_error_2d(matrixv &solution);

    static double exact_solution(double x);

    static double exact_solution(double x,double y);
};

//void printVector(vecd &array) {
//    for (int i = 0; i < array.size(); ++i) {
//        cout << array[i] << '\t' << '\t';
//    }
//    cout << endl;
//}

#endif //MAIN_CPP_ERRORCOMPARE_H
