#include <iostream>
#include <vector>
#include "lib/Function.hpp"

using namespace std;

int main() {
    func::Function f;
    vector<int> t;
    for (int i = 1; i <= 10; i++) {
        t.push_back(i);
    }
    cout << f.mean(t) << endl;
    return 0;
}
