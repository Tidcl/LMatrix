#ifndef MAIN_CPP_MYOBJECT_H
#define MAIN_CPP_MYOBJECT_H

#include "../NMatrixDefine.h"
using namespace std;

class MyPrint {
    int accuracy=0;
public:
    MyPrint(int &accuracy);

    ~MyPrint() = default;

public:

    void print_vector(vecd &vtr) const;

    void print_vector(valarray<double> &valary) const;

    void print_matrix(matrixv &mtx) const;

};


#endif //MAIN_CPP_MYOBJECT_H
