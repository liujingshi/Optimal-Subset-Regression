#ifndef FUNCTION_HPP
#define FUNCTION_HPP
#include <vector>

namespace func
{

class Function
{

  public:
    double MGF(std::vector<double>, int, int);

    std::vector<std::vector<double> > MGF(std::vector<double>, int);

    std::vector<std::vector<double> > ArrIn(std::vector<std::vector<double> >);

    std::vector<double> Differential(std::vector<double>);

    double SumAdd(std::vector<std::vector<double> >, int, int);

    std::vector<std::vector<double> > SumAdd(std::vector<std::vector<double> >, double);

    double CalcQ(std::vector<double>);

    double CalcQ(std::vector<double>, std::vector<double>);

    double CalcS1(std::vector<double>, std::vector<double>);

    double CalcUV(std::vector<double>);

    int CalcT(std::vector<double>, double, int t);

    std::vector<std::vector<double> > CalcNij(std::vector<double>, std::vector<double>);

    double CalcS2(std::vector<double>, std::vector<double>);

    double CalcCSC(std::vector<double>, std::vector<double>);

    double CalcCSC(std::vector<double>, std::vector<std::vector<double> >);

    std::vector<double> CalcCSCs(std::vector<double>, std::vector<std::vector<double> >);

    std::vector<double> CalcCSCs(std::vector<double>, std::vector<std::vector<std::vector<double> > >);

    double CalcX2(std::vector<double>, std::vector<double>);

    std::vector<double> CalcX2(std::vector<double>, std::vector<std::vector<double> >);

    std::vector<int> RouSelect(std::vector<double>, double);

    std::vector<std::vector<double> > RouSelect(std::vector<std::vector<std::vector<double> > >, std::vector<std::vector<double> >, double x);

    void TwoAddOne(std::vector<int>&);

    bool TwoIsFull(std::vector<int>);

    std::vector<std::vector<std::vector<double> > > Group(std::vector<std::vector<double> >);

    std::vector<std::vector<double> > Mul(std::vector<std::vector<double> > , std::vector<std::vector<double> >);

    std::vector<std::vector<double> > Inv(std::vector<std::vector<double> >);

    std::vector<std::vector<double> > Xb(std::vector<double>);

    std::vector<std::vector<double> > Xb(std::vector<std::vector<double> >);

    std::vector<std::vector<double> > ComeBackP(std::vector<double>, std::vector<std::vector<double> >);

    std::vector<double> ComeBack(std::vector<double>, std::vector<std::vector<double> >);

    int FindMaxCSC(std::vector<double>);

    std::vector<double> Predict(std::vector<std::vector<double> >, std::vector<std::vector<double> >, int);

    std::vector<double> Predict(std::vector<double>, int);
};

} // namespace func

#endif