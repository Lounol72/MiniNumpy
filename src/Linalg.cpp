#include "mnp/Linalg.hpp"
#include "mnp/Matrix.hpp"
#include <stdexcept>
#include <vector>
#include <cmath>
#include <utility>

std::vector<double> mnp::solve(const mnp::Matrix &m,
                               const std::vector<double> &v) {
  if (!m.isSquare() || (m.nbLignes() != static_cast<int>(v.size())))
    throw std::invalid_argument("");
  // Copie des paramètres
  mnp::Matrix m_x = m;
  m_x.ajouterColonne(mnp::Vector(v));

  // Pivot de Gauss sur la matrice m_x
  m_x = mnp::pivotDeGauss(m_x);

  // La remontée
  return mnp::remontee(m_x);
}

mnp::Matrix mnp::pivotDeGauss(const mnp::Matrix &m) {
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
      throw std::domain_error(
          "Matrice Singulière : pas de pivot non nul en colonne <" + std::to_string(i) + ">");
    }
    for (int y = i + 1; y < m_x.nbLignes(); y++) {
      double facteur = (m_x(y, i) / pivot);
      for (int j = i; j < m_x.nbColonne(); j++)
        m_x(y, j) = m_x(y, j) - (facteur * m_x(i, j));
    }
  }
  return m_x;
}

std::vector<double> mnp::remontee(const mnp::Matrix &m) {
  std::vector<double> res(m.nbColonne() - 1);

  // Boucle pour l'accès aux valeurs
  for (int i = static_cast<int>(res.size()) - 1; i >= 0; i--) {
    double intermediaire = 0;
    for (int j = i + 1; j < static_cast<int>(res.size()); j++)
      intermediaire -= m(i, j) * res.at(j);

    res.at(i) = (m(i, static_cast<int>(res.size())) +
                 intermediaire) / // derniere colonne
                m(i, i);          // valeur dans la matrice
  }

  return res;
}
