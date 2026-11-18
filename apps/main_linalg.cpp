#include "mnp/Linalg.hpp"

int main(){
  
  const mnp::Matrix systeme({{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}});
  const std::vector<double> vector_equal = {8, -11, 3};
  std::vector<double> v = {2, 3, -1};
  const auto x = mnp::solve(systeme, vector_equal);
  return 0;
}
