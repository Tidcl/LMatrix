//
// Created by Administrater on 2022/12/25.
//

#include <iomanip>
#include "PrintTools.h"

using namespace std;

MyPrint::MyPrint(int &accuracy) {
    this->accuracy = accuracy;
}

void MyPrint::print_vector(vector<double> &vtr) const {
    cout.precision(accuracy);
    cout.setf(ios::fixed);
    string readme = "***----This Vector----***";
    string endme = "**---Show Vector End---**";
    cout << readme << endl;
    cout << "[";
    for (int i = 0; i < vtr.size(); ++i) {
        if (i + 1 == vtr.size()) {
            cout << vtr[i];
        } else
            cout << vtr[i] << ",  ";
    }
    cout << "]" << endl;
    cout << endme << endl;
    cout << "--------------------------------------------------------" << endl;
}

void MyPrint::print_vector(valarray<double> &valary) const {
    cout.precision(accuracy);
    cout.setf(ios::fixed);
    string readme = "***----This Vector----***";
    string endme = "**---Show Vector End---**";
    cout << readme << endl;
    cout << "[";
    for (int i = 0; i < valary.size(); ++i) {
        if (i + 1 == valary.size()) {
            cout << valary[i];
        } else
            cout << valary[i] << ",  ";
    }
    cout << "]" << endl;
    cout << endme << endl;
    cout << "--------------------------------------------------------" << endl;
}

void MyPrint::print_matrix(vector<vector<double>> &mtx) const {
    cout.precision(accuracy);
    cout.setf(ios::fixed);
    string readme = "***----This Matrix----***";
    string endme = "**---Show Matrix End---**";
    cout << readme << endl;
    for (int i = 0; i < mtx.size(); ++i) {
        cout << "[";
        for (int j = 0; j < mtx[0].size(); ++j) {
            if (j + 1 == mtx[0].size()) {
                cout << setw(11) << setiosflags(ios::right) << mtx[i][j];
            } else {
                cout << setw(11) << setiosflags(ios::right) << mtx[i][j] << ",\t";
            }
        }
        cout << " ]" << endl;
    }
    cout << endme << endl;
    cout << "--------------------------------------------------------" << endl;

}


//int main() {
//    vector<double> P = {0, 0.25, 0.5, 0.75, 1};
//    vector<vector<double>> T = {{0, 1, 2, 3},
//                                {1, 2, 3, 4}};
//
//    int accuracy = 6;
//    MyPrint myPrint(accuracy);
//    myPrint.print_matrix(T);
//    myPrint.print_vector(P);
//
//}