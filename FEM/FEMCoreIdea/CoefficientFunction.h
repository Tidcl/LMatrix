//
// Created by admin on 2023/1/5.
//

#ifndef MAIN_CPP_COEFFFUN_H
#define MAIN_CPP_COEFFFUN_H


class CoefficientFun {
private:
    char select_mark;
public:
    CoefficientFun(char select_mark);

    ~CoefficientFun() = default;

    double coefficient_function(double &x, double &c_func);

    double coefficient_function_2d(double &x, double &y,double &c_func);
};


#endif //MAIN_CPP_COEFFFUN_H
