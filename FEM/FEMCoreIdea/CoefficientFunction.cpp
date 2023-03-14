//
// Created by admin on 2023/1/5.
//

#include <cmath>
#include "CoefficientFunction.h"

// 定义了系数函数的函数，传参【选择标记】
CoefficientFun::CoefficientFun(char select_mark) {
    this->select_mark = select_mark;
}

double CoefficientFun::coefficient_function(double &x, double &c_func) {
    if (select_mark == 'c') {  // 对应刚度矩阵部分
        c_func = exp(x);
        return c_func;
    } else if (select_mark == 'f') {  // 对应载荷向量部分
        c_func = -exp(x) * (cos(x) - 2 * sin(x) - x * cos(x) - x * sin(x));
        return c_func;
    } else if (select_mark == 'b') {  // 留着改善测试用
        c_func = cos(x);
        return c_func;
    } else {
        return 0;
    }
}

double CoefficientFun::coefficient_function_2d(double &x, double &y, double &c_func) {
    if (select_mark == 'c') {  // 对应刚度矩阵部分
        c_func = pow(x,0);
        return c_func;
    } else if (select_mark == 'f') {  // 对应载荷向量部分
        c_func = -y * (1 - y) * (1 - x - pow(x , 2) / 2) * exp(x + y) - x * (1 - x / 2) * (
                -3 * y - pow(y , 2)) * exp(x + y);
        return c_func;
    } else if (select_mark == 'b') {  // 留着改善测试用
        c_func = cos(x);
        return c_func;
    } else {
        return 0;
    }
}