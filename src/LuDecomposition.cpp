#include "mnp/LuDecomposition.hpp"
#include <stdexcept>

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
} // namespace mnp
