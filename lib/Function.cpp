#include <vector>
#include <stdio.h>
#include "Function.hpp"

double func::Function::MGF(std::vector<double> x, int l, int i) {
    double sum = 0;
    int num = (int)(x.size() / l);
    for (int k = i - 1; k < x.size(); k += l) {
        sum += x[k];
    }
    return sum / num;
}

std::vector<std::vector<double> > func::Function::MGF(std::vector<double> x) {
    std::vector<std::vector<double> > F;
    for (int l = 1; l <= x.size(); l++) {
        std::vector<double> H;
        while (H.size() < x.size()) {
            for (int i = 1; i <= l && H.size() < x.size(); i++) {
                H.push_back(this->MGF(x, l, i));
            }
            
        }
        F.push_back(H);
    }
    return this->ArrIn(F);
}

std::vector<std::vector<double> > func::Function::ArrIn(std::vector<std::vector<double> > F) {
    for (int i = 0; i < F.size(); i++) {
        for (int j = 0; j < F[i].size(); j++) {
            if (i > j) {
                double tem = F[i][j];
                F[i][j] = F[j][i];
                F[j][i] = tem;
            }
        }
    }
    return F;
}

std::vector<double> func::Function::Differential(std::vector<double> x) {
    std::vector<double> xd;
    for (int i = 1; i < x.size(); i++) {
        xd.push_back(x[i] - x[i-1]);
    }
    return xd;
}

std::vector<std::vector<double> > func::Function::SumAdd(std::vector<std::vector<double> > x) {
    std::vector<std::vector<double> > xa(x);
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x.size(); j++) {
            
        }
    }
    return xa;
}

void func::Function::show(std::vector<std::vector<double> > F) {
    for (int i = 0; i < F.size(); i++) {
        for (int j = 0; j < F[i].size(); j++) {
            printf("%5.1lf", F[i][j]);
        }
        printf("\n");
    }
}
