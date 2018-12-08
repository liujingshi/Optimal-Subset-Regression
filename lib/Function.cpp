#include <vector>
#include <stdio.h>
#include <cmath>
#include "Function.hpp"

double func::Function::MGF(std::vector<double> x, int l, int i)
{
    double sum = 0;
    int num = (int)(x.size() / l);
    for (int k = i - 1; k < x.size(); k += l)
    {
        sum += x[k];
    }
    return sum / num;
}

std::vector<std::vector<double> > func::Function::MGF(std::vector<double> x)
{
    std::vector<std::vector<double> > F;
    for (int l = 1; l <= x.size(); l++)
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
    return this->ArrIn(F);
}

std::vector<std::vector<double> > func::Function::ArrIn(std::vector<std::vector<double> > F)
{
    for (int i = 0; i < F.size(); i++)
    {
        for (int j = 0; j < F[i].size(); j++)
        {
            if (i > j)
            {
                double tem = F[i][j];
                F[i][j] = F[j][i];
                F[j][i] = tem;
            }
        }
    }
    return F;
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
        for (int i = 1; i < xa.size(); i++)
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

