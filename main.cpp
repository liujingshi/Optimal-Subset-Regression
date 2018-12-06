
#include <iostream>
#include "lib/Hello.hpp"

using namespace std;

int main() {
    Hello *h = new Hello;
    cout << h->get() << endl;
    return 0;
}
