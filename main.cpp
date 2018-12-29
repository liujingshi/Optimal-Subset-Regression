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
    double a[100] = {7, 4, 0, 0, 4, 9, 2, 3, 4, 8, 1, 1, 8, 9, 9, 6, 5, 7, 7, 4, 1, 9, 3, 5, 9, 7, 7, 8, 4, 2, 1, 7, 4, 4, 6, 5, 2, 5, 7, 3, 1, 9, 1, 1, 7, 6, 3, 0, 0, 5, 6, 7, 2, 0, 7, 1, 6, 9, 7, 8, 8, 3, 0, 6, 5, 7, 0, 4, 9, 1, 4, 2, 1, 5, 1, 8, 8, 8, 9, 7, 6, 1, 8, 4, 8, 9, 7, 2, 7, 2, 6, 3, 9, 1, 2, 6, 2, 9, 9, 7};
    vector<double> x(a, a + 100);
    //vector<double> result = f.Predict(x, 1);
    //show(result);
    vector<double> result = f1.PredictE(x, 1);
    show(result);
    return 0;
}
