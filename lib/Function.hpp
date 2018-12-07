#ifndef FUNCTION_HPP
#define FUNCTION_HPP
#include <vector>

namespace func {

    class Function {

    public:

        double MGF(std::vector<double>, int, int);

        std::vector<std::vector<double> > MGF(std::vector<double>);

        std::vector<std::vector<double> > ArrIn(std::vector<std::vector<double> >);

        std::vector<double> Differential(std::vector<double>);

        double SumAdd(std::vector<std::vector<double> >, int, int);

        std::vector<std::vector<double> > SumAdd(std::vector<std::vector<double> >, double);

        double CalcQ(std::vector<double>);

        double CalcQ(std::vector<double>, std::vector<double>);

        double CalcS1(std::vector<double>, std::vector<double>);

        //double CalcS2();

        //double CalcCSC();

        void show(std::vector<std::vector<double> >);

    };

}

#endif