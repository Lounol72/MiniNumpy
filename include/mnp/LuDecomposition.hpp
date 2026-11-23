#pragma once

#include "mnp/Matrix.hpp"
#include <vector>
namespace mnp {

class LuDecomposition {
public:
  // Constructeur
  LuDecomposition(const mnp::Matrix m)
      : m_mat(m.nbLignes(), m.nbColonne()), m_p(m.nbLignes(), 0) {
    mnp::Matrix m_x = m;

    std::for_each(m_p.begin(), m_p.end(),
                  [idx = 0](int &elem) mutable { elem = idx++; });

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

      // Echange si besoin
      if (lignePivot != i) {
        // swap lignes
        m_nbEchanges++;
        for (int y = 0; y < m_x.nbColonne(); y++) {
          std::swap(m_x(i, y), m_x(lignePivot, y));
        }
        std::swap(m_p[i], m_p[lignePivot]);
      }

      // Vérification de matrice singulière
      if (approx_equal(pivot, 0)) {
        m_singuliere = true;
        continue;
      }
      // Application du pivot avec le facteur
      for (int y = i + 1; y < m_x.nbLignes(); y++) {
        double facteur = (m_x(y, i) / pivot);
        for (int j = i; j < m_x.nbColonne(); j++)
          m_x(y, j) = m_x(y, j) - (facteur * m_x(i, j));
        // sauvegarde du facteur
        m_x(y, i) = facteur;
      }
    }
    m_mat = m_x;
  }

  // Copie 
  LuDecomposition(const LuDecomposition &) = default;
  LuDecomposition(LuDecomposition &&) = default;
  LuDecomposition &operator=(const LuDecomposition &) = delete;
  LuDecomposition &operator=(LuDecomposition &&) = delete;

  // Accesseur
  double determinant() const;
  std::vector<double> resoudre(std::vector<double> b) const;
  mnp::Matrix inverse() const;

private:
  mnp::Matrix m_mat;
  std::vector<int> m_p;
  bool m_singuliere = false;
  int m_nbEchanges = 0;
};
} // namespace mnp
