#ifndef MAIN_CPP_MYOBJECT_H
#define MAIN_CPP_MYOBJECT_H

#include <iostream>
#include <vector>
#include <valarray>

using namespace std;

class MyPrint {
    int accuracy=0;
public:
    MyPrint(int &accuracy);

    ~MyPrint() = default;

public:

    void print_vector(vector<double> &vtr) const;

    void print_vector(valarray<double> &valary) const;

    void print_matrix(vector<vector<double>> &mtx) const;

};


#endif //MAIN_CPP_MYOBJECT_H
