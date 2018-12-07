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

        std::vector<std::vector<double> > SumAdd(std::vector<std::vector<double> >);

        void show(std::vector<std::vector<double> >);

    };

}

#endif