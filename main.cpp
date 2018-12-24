#include <iostream>
#include <vector>
#include <stdio.h>
#include "lib/Function.hpp"

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
    func::Function f;
    double a[36] = {38.3,35.9,36.2,37.5,37.0,40.5,39.0,36.6,36.6,36.6,38.1,36.7,35.7,36.0,36.0,37.6,36.7,39.7,36.5,36.0,37.9,34.4,38.2,36.4,36.8,36.7,35.2,38.5,36.4,36.8,37.0,37.9,35.7,38.0,37.2,36.5};
    vector<double> x(a, a + 36);
    vector<double> result = f.Predict(x, 5);
    show(result);
    return 0;
}
