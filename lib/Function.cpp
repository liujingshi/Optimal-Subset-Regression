#include <vector>
#include "Function.hpp"

int func::Function::mean(std::vector<int> t) {
    int sum = 0;
    int num = t.size();
    while (!t.empty()) {
        sum += t.back();
        t.pop_back();
    }
    return (int)(sum / num);
}
