//
// Created by admin on 2023/1/5.
//
#include <iostream>
#include<vector>

using namespace std;

//vector<vector<double>> P = {{0, 0.25, 0.5, 0.75, 1}};
//vector<vector<double>> T = {{0, 1, 2, 3},
//                         {1, 2, 3, 4}};

vector<vector<double>> P = {{0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0}};
vector<vector<double>> T = {{0, 2, 4, 6},
                            {1, 3, 5, 7},
                            {2, 4, 6, 8}};

void print_matrix_asm(vector<vector<double>> &matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[0].size(); ++j) {
            cout << matrix[i][j] << '\t' << '\t';
        }
        cout << endl;
    }
}

int main() {
    cout << P[0][int(T[0][0])] << endl;
    vector<vector<double>> vertices(1, vector<double>(T.size(), 0.0));
    print_matrix_asm(vertices);
    cout << "----" << endl;
    for (int n = 0; n < 4; ++n) {
//        vertices = P[:,T[:,n]];
        for (int i = 0; i < P.size(); ++i) {
            for (int j = 0; j < T.size(); ++j) {
                for (int k = 0; k < n + 1; ++k) {
                    vertices[i][j] = P[i][int(T[j][k])];
                }
            }
        }
        print_matrix_asm(vertices);
    }
}
