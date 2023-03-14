#include <iostream>
#include <vector>
#include "LinearEqn.h"

using namespace std;

int main() {
    cout.setf(ios::fixed);
    cout.precision(4);

    vector<double> b = {11, 9, 9, 9, 13, 17};
    vector<vector<double>> a = {{6, 5, 4, 3, 2, 1},
                                {5, 6, 5, 4, 3, 2},
                                {4, 5, 6, 5, 4, 3},
                                {3, 4, 5, 6, 5, 4},
                                {2, 3, 4, 5, 6, 5},
                                {1, 2, 3, 4, 5, 6}};
//    solution = {3,-1,0,-2,0,4};

    double tol = 1e-14;
    int maxIter = 1000;
    int m = 6;
    int n = 6;

    LinearEqn L1(m, n, a, b);
    L1.BicgStab(tol, maxIter);
    L1.printSolution();
//    L1.printVector(b);

    return 0;
}