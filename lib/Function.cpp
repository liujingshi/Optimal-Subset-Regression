#include <vector>
#include <cmath>
#include <stdio.h>
#include "Function.hpp"

double func::Function::MGF(std::vector<double> x, int l, int i)
{
    double sum = 0;
    int num = (int)(x.size() / l);
    for (int j = 0; j < num; j++) {
        sum += x[i - 1 + (j * l)];
    }
    return sum / num;
}

std::vector<std::vector<double> > func::Function::MGF(std::vector<double> x)
{
    std::vector<std::vector<double> > F;
    int m = (int)(x.size() / 3);
    for (int l = 1; l <= m; l++)
    {
        std::vector<double> H;
        while (H.size() < x.size())
        {
            for (int i = 1; i <= l && H.size() < x.size(); i++)
            {
                H.push_back(this->MGF(x, l, i));
            }
        }
        F.push_back(H);
    }
    return F;
}

std::vector<std::vector<double> > func::Function::ArrIn(std::vector<std::vector<double> > F)
{
    int w = F.size();
    int h = F[0].size();
    std::vector<double> item(w, 0);
    std::vector<std::vector<double> > Arr(h, item);
    for (int i = 0; i < F.size(); i++)
    {
        for (int j = 0; j < F[i].size(); j++)
        {
            Arr[j][i] = F[i][j];
        }
    }
    return Arr;
}

std::vector<double> func::Function::Differential(std::vector<double> x)
{
    std::vector<double> xd;
    for (int i = 1; i < x.size(); i++)
    {
        xd.push_back(x[i] - x[i - 1]);
    }
    return xd;
}

double func::Function::SumAdd(std::vector<std::vector<double> > x, int l, int t)
{
    double sum = 0;
    for (int i = 0; i < t; i++)
    {
        sum += x[l][i + 1];
    }
    return sum;
}

std::vector<std::vector<double> > func::Function::SumAdd(std::vector<std::vector<double> > x, double x1)
{
    std::vector<std::vector<double> > xa(x);
    for (int i = 0; i < xa.size(); i++)
    {
        xa[i][0] = x1;
    }
    for (int l = 0; l < xa.size(); l++)
    {
        for (int i = 1; i < xa[l].size(); i++)
        {
            xa[l][i] = x1 + this->SumAdd(x, l, i);
        }
    }
    return xa;
}

double func::Function::CalcQ(std::vector<double> x)
{
    double xb = this->MGF(x, 1, 1);
    int n = x.size();
    double sum = 0;
    for (int t = 0; t < n; t++)
    {
        sum += (x[t] - xb) * (x[t] - xb);
    }
    return sum / n;
}

double func::Function::CalcQ(std::vector<double> x, std::vector<double> xt)
{
    int n = xt.size();
    double sum = 0;
    for (int t = 0; t < n; t++)
    {
        sum += (x[t] - xt[t]) * (x[t] - xt[t]);
    }
    return sum / n;
}

double func::Function::CalcS1(std::vector<double> x, std::vector<double> xk)
{
    double Qx = this->CalcQ(x);
    double Qk = this->CalcQ(x, xk);
    int n = x.size();
    double S1 = n * (1 - (Qk / Qx));
    return S1;
}

double func::Function::CalcUV(std::vector<double> x)
{
    double sum = 0;
    double n = x.size();
    for (int i = 0; i < 0; i++)
    {
        sum += fabs(x[i]);
    }
    return sum / n;
}

int func::Function::CalcT(std::vector<double> x, double u, int t)
{
    int T = 0;
    if (x[t] > u)
    {
        T = 0;
    }
    if (fabs(x[t]) <= u)
    {
        T = 1;
    }
    if (x[t] < -u)
    {
        T = 2;
    }
    return T;
}

std::vector<std::vector<double> > func::Function::CalcNij(std::vector<double> x, std::vector<double> f)
{
    std::vector<double> x1 = this->Differential(x);
    std::vector<double> f1 = this->Differential(f);
    double U = this->CalcUV(x1);
    double V = this->CalcUV(f1);
    std::vector<double> Nj(3, 0);
    int i = 0, j = 0;
    std::vector<std::vector<double> > Nij(3, Nj);
    for (int t = 0; t < f1.size(); t++) {
        i = this->CalcT(x1, U, t);
        j = this->CalcT(f1, V, t);
        Nij[i][j]++;
    }
    return Nij;
}

double func::Function::CalcS2(std::vector<double> x, std::vector<double> f) {
    std::vector<std::vector<double> > Nij = this->CalcNij(x, f);
    double R1 = 0, R2 = 0, R3 = 0, Ni = 0, Nj = 0, tem = 0;
    int n = x.size();
    for (int i = 0; i < Nij.size(); i++) {
        tem = 0;
        for (int j = 0; j < Nij[i].size(); j++) {
            if (Nij[i][j] > 0) {
                tem += Nij[i][j] * log(Nij[i][j]);
            }
        }
        R1 += tem;
    }
    for (int i = 0; i < Nij.size(); i++) {
        Ni = 0;
        for (int j = 0; j < Nij.size(); j++) {
            Ni += Nij[i][j];
        }
        if (Ni > 0) {
            R2 += Ni * log(Ni);
        }
    }
    for (int j = 0; j < Nij.size(); j++) {
        Nj = 0;
        for (int i = 0; i < Nij.size(); i++) {
            Nj += Nij[i][j];
        }
        if (Nj > 0) {
            R3 += Nj * log(Nj);
        }
    }
    return 2 * (R1 + (n - 1) * log(n - 1) - R2 - R3);
}

double func::Function::CalcCSC(std::vector<double> x, std::vector<double> f) {
    double S1 = this->CalcS1(x, f);
    double S2 = this->CalcS2(x, f);
    return S1 + S2;
}

std::vector<double> func::Function::CalcCSC(std::vector<double> x, std::vector<std::vector<double> > F) {
    std::vector<double> CSC;
    for (int i = 0; i < F.size(); i++) {
        CSC.push_back(this->CalcCSC(x, F[i]));
    }
    return CSC;
}

double func::Function::CalcX2(std::vector<double> x, std::vector<double> f) {
    double X2 = 0;
    for (int i = 0; i < f.size(); i++) {
        if (f[i] != 0) {
            X2 += (x[i] - f[i]) * (x[i] - f[i]) / f[i];
        }
    }
    return X2;
}

std::vector<double> func::Function::CalcX2(std::vector<double> x, std::vector<std::vector<double> > F) {
    std::vector<double> X2;
    for (int i = 0; i < F.size(); i++) {
        X2.push_back(this->CalcX2(x, F[i]));
    }
    return X2;
}

std::vector<int> func::Function::RouSelect(std::vector<double> c, double x) {
    std::vector<int> div;
    for (int i = 0; i < c.size(); i++) {
        if (c[i] > x) {
            div.push_back(i);
        }
    }
    return div;
}

std::vector<std::vector<double> > func::Function::RouSelect(std::vector<std::vector<std::vector<double> > > F, std::vector<std::vector<double> > CSC, double x) {
    std::vector<std::vector<double> > P;
    for (int i = 0; i < CSC.size(); i++) {
        std::vector<int> temCSC = this->RouSelect(CSC[i], x);
        for (int j = 0; j < temCSC.size(); j++) {
            P.push_back(F[i][temCSC[j]]);
        }
    }
    return P;
}

std::vector<std::vector<double> > func::Function::Group(std::vector<std::vector<double> > P) {
    std::vector<std::vector<double> > Son;
    
}
