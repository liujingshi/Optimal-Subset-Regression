#include <iostream>
#include <vector>
#include <stdio.h>
#include "lib/Function.hpp"

using namespace std;

void show(vector<vector<double> > F)
{
    for (int i = 0; i < F.size(); i++)
    {
        for (int j = 0; j < F[i].size(); j++)
        {
            printf("%5.1lf", F[i][j]);
        }
        printf("\n");
    }
}

void show(vector<double> x)
{
    for (int i = 0; i < x.size(); i++)
    {
        printf("%8.2lf", x[i]);
    }
    printf("\n");
}

int main()
{
    func::Function f;
    double a[36] = {38.3,35.9,36.2,37.5,37.0,40.5,39.0,36.6,36.6,36.6,38.1,36.7,35.7,36.0,36.0,37.6,36.7,39.7,36.5,36.0,37.9,34.4,38.2,36.4,36.8,36.7,35.2,38.5,36.4,36.8,37.0,37.9,35.7,38.0,37.2,36.5};
    vector<double> x(a, a + 36);
    vector<vector<double> > F0 = f.MGF(x);
    vector<double> x1 = f.Differential(x);
    vector<vector<double> > F1 = f.MGF(x1);
    vector<double> x2 = f.Differential(x1);
    vector<vector<double> > F2 = f.MGF(x2);
    vector<vector<double> > F3 = f.SumAdd(F1, x[0]);
    show(F0);
    show(F1);
    show(F2);
    show(F3);
    vector<double> CSC0 = f.CSC(x, F0);
    vector<double> CSC1 = f.CSC(x, F1);
    vector<double> CSC2 = f.CSC(x, F2);
    vector<double> CSC3 = f.CSC(x, F3);
    /*vector<double> X20 = f.CalcX2(x, F0);
    vector<double> X21 = f.CalcX2(x, F1);
    vector<double> X22 = f.CalcX2(x, F2);
    vector<double> X23 = f.CalcX2(x, F3);*/
    show(CSC0);
    show(CSC1);
    show(CSC2);
    show(CSC3);
    /*show(X20);
    show(X21);
    show(X22);
    show(X23);*/
    return 0;
}
