#include "mnp/LuDecomposition.hpp"
#include <stdexcept>
#include <vector>

namespace mnp {
std::vector<double>
LuDecomposition::resoudre(const std::vector<double> b) const {
  if (m_singuliere)
    throw std::domain_error(
        "Résolution impossible pour une matrice singuliere");

  if (static_cast<int>(b.size()) != m_mat.nbLignes())
    throw std::invalid_argument("b n'est pas de la même taille que A");

  std::vector<double> b_x(b.size());
  // Pb : on applique à b les échanges de lignes faits sur A (PA = LU)
  for (size_t i = 0; i < m_p.size(); i++) {
    b_x[i] = b[m_p[i]];
  }
  std::vector<double> res_x(m_mat.nbLignes());
  // Descente : L y = Pb (diagonale de L implicite à 1, donc pas de division)
  for (int i = 0; i < static_cast<int>(res_x.size()); i++) {
    double intermediaire = 0;
    for (int j = 0; j < i; j++) {
      intermediaire -= m_mat(i, j) * res_x.at(j);
    }
    res_x.at(i) = (b_x.at(i) + intermediaire);
  }

  // Remontée : U x = y (diagonale de U stockée dans m_mat)
  std::vector<double> res(m_mat.nbLignes());
  for (int i = static_cast<int>(res.size()) - 1; i >= 0; i--) {
    double intermediaire = 0;
    for (int j = i + 1; j < static_cast<int>(res.size()); j++) {
      intermediaire -= m_mat(i, j) * res.at(j);
    }

    res.at(i) = (res_x.at(i) + intermediaire) / m_mat(i, i);
  }

  return res;
}

mnp::Matrix LuDecomposition::resoudreMatrice(const mnp::Matrix &b) const {
  if (m_singuliere)
    throw std::domain_error(
        "Résolution impossible pour une matrice singuliere");

  if (static_cast<int>(b.nbLignes()) != m_mat.nbLignes())
    throw std::invalid_argument("b n'as pas le même nombre de lignes que A");

  mnp::Matrix b_x(b.nbLignes(), b.nbColonne());

  // Pb : on applique à b les échanges de lignes faits sur A (PA = LU)
  for (int j = 0; j < b.nbColonne(); j++) {
    for (size_t i = 0; i < m_p.size(); i++) {
      b_x(i, j) = b(m_p[i], j);
    }
  }

  // Descente : L Y = PB, sur place dans b_x (diagonale de L implicite à 1)
  for (int i = 0; i < m_mat.nbLignes(); i++) {
    for (int j = 0; j < i; j++) {
      const double l = m_mat(i, j); // lu une seule fois
      for (int k = 0; k < b.nbColonne(); k++)
        b_x(i, k) -= l * b_x(j, k); // appliqué aux K seconds membres
    }
  }

  // Remontée : U X = Y, sur place dans b_x (diagonale de U stockée dans m_mat)
  for (int i = m_mat.nbLignes() - 1; i >= 0; i--) {
    for (int j = i + 1; j < m_mat.nbColonne(); j++) {
      const double u = m_mat(i, j); // lu une seule fois
      for (int k = 0; k < b.nbColonne(); k++)
        b_x(i, k) -= u * b_x(j, k); // lignes j > i : déjà résolues
    }
    // la ligne i est complète : on divise par le pivot une fois pour toutes
    const double pivot = m_mat(i, i);
    for (int k = 0; k < b.nbColonne(); k++)
      b_x(i, k) /= pivot;
  }

  return b_x;
}
} // namespace mnp
