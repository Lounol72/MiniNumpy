#include <iostream>

#include "mnp/version.hpp"

int main() {
    std::cout << "mini-numpy " << mnp::version() << '\n';
    return 0;
}
