#pragma once

#include "mnp/Matrix.hpp"
#include "mnp/Precision.hpp"
#include <stdexcept>
#include <vector>
namespace mnp {
inline const mnp::Matrix &carree(const mnp::Matrix &m) {
  if (!m.isSquare())
    throw std::invalid_argument("LuDecomposition : matrice non carrée");
  return m;
}
class LuDecomposition {
public:
  // Constructeur
  explicit LuDecomposition(const mnp::Matrix &m)
      : m_mat(carree(m)), m_p(m.nbLignes(), 0) {
    std::for_each(m_p.begin(), m_p.end(),
                  [idx = 0](int &elem) mutable { elem = idx++; });

    for (int i = 0; i < m_mat.nbLignes(); i++) {
      int lignePivot = i;
      double pivot = m_mat(i, i);

      // Trouver le plus grand pivot
      for (int y = i + 1; y < m_mat.nbLignes(); y++) {
        if (std::fabs(m_mat(y, i)) > std::fabs(pivot)) {
          lignePivot = y;
          pivot = m_mat(y, i);
        }
      }

      // Echange si besoin
      if (lignePivot != i) {
        // swap lignes
        m_nbEchanges++;
        for (int y = 0; y < m_mat.nbColonne(); y++) {
          std::swap(m_mat(i, y), m_mat(lignePivot, y));
        }
        std::swap(m_p[i], m_p[lignePivot]);
      }

      // Vérification de matrice singulière
      if (approx_equal(pivot, 0)) {
        m_singuliere = true;
        continue;
      }
      // Application du pivot avec le facteur
      for (int y = i + 1; y < m_mat.nbLignes(); y++) {
        double facteur = (m_mat(y, i) / pivot);
        for (int j = i; j < m_mat.nbColonne(); j++)
          m_mat(y, j) = m_mat(y, j) - (facteur * m_mat(i, j));
        // sauvegarde du facteur
        m_mat(y, i) = facteur;
      }
    }
  }

  // Copie
  LuDecomposition(const LuDecomposition &) = default;
  LuDecomposition(LuDecomposition &&) = default;
  LuDecomposition &operator=(const LuDecomposition &) = delete;
  LuDecomposition &operator=(LuDecomposition &&) = delete;

  // Accesseur
  double determinant() const;
  std::vector<double> resoudre(const std::vector<double> b) const;
  mnp::Matrix inverse() const;
  bool estSinguliere() const { return m_singuliere; }

private:
  mnp::Matrix m_mat;
  std::vector<int> m_p;
  bool m_singuliere = false;
  int m_nbEchanges = 0;
};
} // namespace mnp
