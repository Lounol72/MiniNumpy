#include "mnp/Linalg.hpp"
#include "mnp/Matrix.hpp"
#include <iostream>
#include <vector>

std::vector<double> mnp::solve(const mnp::Matrix &m,
                               const std::vector<double> &v) {
  // Copie des paramètres
  mnp::Matrix m_x(m.nbLignes(), m.nbColonne() + 1);
  // On rempli la matrice augmentée
  for (int i = 0; i < m.nbLignes(); i++) {
    for (int j = 0; j < m.nbColonne(); j++)
      m_x(i, j) = m(i, j);
    m_x(i, m_x.nbColonne() - 1) = v[i];
  }

  for (int i = 0; i < static_cast<int>(v.size()); i++) {
    double pivot = m_x(i, i);
    for (int y = i + 1; y < static_cast<int>(v.size()); y++) {
      double facteur = (m_x(y, i) / pivot);
      for (int j = 0; j < m_x.nbColonne(); j++)
        m_x(y, j) = m_x(y, j) - (facteur * m_x(i, j));
    }
  }

  // La remontée
  return mnp::remontee(m_x);
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
