#pragma once

#include "mnp/Matrix.hpp"
namespace mnp {

class LuDecomposition {
public:
LuDecomposition(const mnp::Matrix m) : mat(m.nbLignes(), m.nbColonne()) {
    mnp::Matrix m_x = m;
    for (int i = 0; i < m_x.nbLignes(); i++) {
      int lignePivot = i;
      double pivot = m_x(i, i);
      // Trouver le plus grand pivot
      for (int y = i + 1; y < m_x.nbLignes(); y++) {
        if (std::fabs(m_x(y, i)) > std::fabs(pivot)) {
          lignePivot = y;
          pivot = m_x(y, i);
        }
      }
      if (lignePivot != i) {
        // swap lignes
        for (int y = 0; y < m_x.nbColonne(); y++)
          std::swap(m_x(i, y), m_x(lignePivot, y));
      }
      // Vérification de matrice singulière
      if (approx_equal(pivot, 0)) {
        singuliere = true;
      }
      for (int y = i + 1; y < m_x.nbLignes(); y++) {
        double facteur = (m_x(y, i) / pivot);
        mat(y,i) = facteur;
        for (int j = i; j < m_x.nbColonne(); j++)
          m_x(y, j) = m_x(y, j) - (facteur * m_x(i, j));
      }
    }
    mat +=m_x;
  }

  LuDecomposition(const LuDecomposition &) = default;
  LuDecomposition(LuDecomposition &&) = default;
  LuDecomposition &operator=(const LuDecomposition &) = delete;
  LuDecomposition &operator=(LuDecomposition &&) = delete;

private:
  mnp::Matrix mat;
  bool singuliere = false;
};
} // namespace mnp
