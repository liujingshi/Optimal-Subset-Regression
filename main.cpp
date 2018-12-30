#include <iostream>
#include <vector>
#include <stdio.h>
//#include "lib/Function.hpp"
#include "lib/Function.1.hpp"

using namespace std;

void show(vector<vector<vector<double> > > F);
void show(vector<vector<double> > F);
void show(vector<double> x);

void show(vector<vector<vector<double> > > F) {
    printf("[");
    for (int i = 0; i < F.size(); i++)
    {
        show(F[i]);
        if (i < F.size() - 1) {
            printf(", ");
        }
        printf("\n");
    }
    printf("]\n");
}

void show(vector<vector<double> > F)
{
    printf("[");
    for (int i = 0; i < F.size(); i++)
    {
        show(F[i]);
        if (i < F.size() - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

void show(vector<double> x)
{
    printf("[");
    for (int i = 0; i < x.size(); i++)
    {
        printf("%.1lf", x[i]);
        if (i < x.size() - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

int main()
{
    //func::Function f;
    func1::Function f1;
    double a[22] = {10, 11, 9, 12, 10, 13, 11, 14, 15, 14, 16, 15, 17, 18, 19, 21, 19, 20, 18, 21, 19, 22};
    vector<double> x(a, a + 22);
    //vector<double> result = f.Predict(x, 1);
    //show(result);
    vector<double> result = f1.PredictE(x, 10);
    show(result);
    return 0;
}
