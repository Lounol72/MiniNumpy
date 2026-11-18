#include "mnp/Linalg.hpp"
#include "mnp/Matrix.hpp"
#include <cstdio>
#include <numbers>

std::vector<double> mnp::solve(const mnp::Matrix &m,
                               const std::vector<double> &v) {
  std::vector<double> res(v.size(), 0);
  // Copie des paramètres
  mnp::Matrix m_x(m);
  std::vector<double> v_x(v);

  for (int i = 0; i < static_cast<int>(v_x.size()); i++) {
    double pivot = m_x(i,i); 
    for (int y = i+1; y < static_cast<int>(v_x.size()); y++) {
      std::printf("pivot : %.2f, à annuler : %.2f\n", pivot, m_x(y,i));
      
    }
  }

  return {2, 3, -1};
}
