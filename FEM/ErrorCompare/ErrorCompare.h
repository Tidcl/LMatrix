//
// Created by admin on 2023/1/7.
//

#ifndef MAIN_CPP_ERRORCOMPARE_H
#define MAIN_CPP_ERRORCOMPARE_H

#include <iostream>
#include <vector>

using namespace std;

class ErrorCompareAnalysis {
private:
    vector<vector<double>> P;
    int num_node;
//    int left;
//    double h_basis;

public:
    ErrorCompareAnalysis(int &num_node, vector<vector<double>> &P);

    ~ErrorCompareAnalysis() = default;

public:

    double get_entire_error_1d(vector<vector<double>> &solution,
                               vector<double> &error_list, bool show_error_list) const;

    double get_entire_error_2d(vector<vector<double>> &solution);

    static double exact_solution(double x);

    static double exact_solution(double x,double y);
};

//void printVector(vector<double> &array) {
//    for (int i = 0; i < array.size(); ++i) {
//        cout << array[i] << '\t' << '\t';
//    }
//    cout << endl;
//}

#endif //MAIN_CPP_ERRORCOMPARE_H
