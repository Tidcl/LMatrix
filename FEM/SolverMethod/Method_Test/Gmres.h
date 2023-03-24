#ifndef GAUSSEPP_CPP_GMRES_H
#define GAUSSEPP_CPP_GMRES_H

#include <iostream>
#include "../../NMatrixDefine.h"

using namespace std;

template<class Tv = vecd>
class Gmres {
public:
    int m_size{};
    Tv m_b;
    vector<Tv> m_mtx;
public:
    Gmres();
    Gmres(const Gmres &) = default;
    Gmres(int);
    Gmres(int, int, const vector<Tv> &m_mtx, const Tv &m_b);

    ~Gmres();
public:
    void printMatrix(vector<Tv> &);

    void printVector(Tv &);

    bool solve(int);

    void Arnoldi(int m, vector<Tv> &h, vector<Tv> &v, double beta);

    void Givens(int m, vector<Tv> &h, Tv &g, double beta,Tv &x, vector<Tv> &v);

public:
    friend double dot(const Tv &, const Tv &);

    friend double norm2(const Tv &);

    friend Tv operator*(double C, const Tv &);

    friend Tv operator*(const vector<Tv> &, const Tv &);

    friend Tv operator-(const Tv &, const Tv &);
};


#endif //GAUSSEPP_CPP_GMRES_H
