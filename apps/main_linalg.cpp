#include "mnp/Linalg.hpp"
#include <iostream>

int main() {

  const mnp::Matrix systeme({{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}});
  const std::vector<double> vector_equal = {8, -11, -3};
  std::vector<double> v = {2, 3, -1};
  const auto x = mnp::solve(systeme, vector_equal);
  for(const auto &val:x)
    std::cout<< val << ", ";
  std::cout << std::endl;
  return 0;
}
