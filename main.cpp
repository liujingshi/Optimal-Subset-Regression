#include <iostream>
#include <vector>
#include <stdio.h>
#include "lib/Function.hpp"

using namespace std;

int main() {
    func::Function f;
    double a[21] = {27.6, 26.8, 27.7, 28.0, 28.0, 27.4, 26.8, 26.9, 28.1, 28.0, 28.0, 27.8, 27.9, 27.2, 27.7, 26.8, 28.0, 27.6, 27.3, 26.9, 28.1};
    vector<double> x(a, a+21);
    vector<vector<double> > F0 = f.MGF(x);
    vector<double> x1 = f.Differential(x);
    vector<vector<double> > F1 = f.MGF(x1);
    vector<double> x2 = f.Differential(x1);
    vector<vector<double> > F2 = f.MGF(x2);
    vector<vector<double> > F3 = f.SumAdd(F1, x[0]);
    f.show(F0);
    f.show(F1);
    f.show(F2);
    f.show(F3);
    return 0;
}
