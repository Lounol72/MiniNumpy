#include "mnp/Matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mnp {

// ==========================================================================
// Constructeurs
// ==========================================================================

Matrix::Matrix(int i, int j) : ligne(i), colonne(j), m_values(i * j, 0.0) {}

Matrix::Matrix(const std::vector<std::vector<double>> &vect)
    : ligne(vect.size()), colonne(vect.empty() ? 0 : vect[0].size()),
      m_values(colonne * ligne, 0.0) {
  if (ligne == 0 || colonne == 0)
    throw std::invalid_argument("Error: empty matrix");

  for (int i = 0; i < ligne; i++) {
    if (static_cast<int>(vect[i].size()) != colonne)
      throw std::invalid_argument("Error: all lines must have the same length");
    for (int j = 0; j < colonne; j++)
      m_values[idx(i, j)] = vect[i][j];
  }
}

Matrix Matrix::identity(int n) {
  Matrix m(n, n);
  for (int i = 0; i < n; i++)
    m.m_values[m.idx(i, i)] = 1.0;
  return m;
}

// ==========================================================================
// Propriétés
// ==========================================================================

bool Matrix::isSquare() const { return ligne == colonne; }

bool Matrix::isZero() const {
  return std::all_of(m_values.begin(), m_values.end(),
                     [](double x) { return std::fabs(x) < PRECISION; });
}

bool Matrix::isIndentity() const {
  if (!isSquare())
    return false;

  for (int i = 0; i < ligne; i++)
    for (int j = 0; j < colonne; j++) {
      const double attendu = (i == j) ? 1.0 : 0.0;
      if (std::fabs(m_values[idx(i, j)] - attendu) >= PRECISION)
        return false;
    }
  return true;
}

Matrix Matrix::transpose() const {
  Matrix t(colonne, ligne);
  for (int i = 0; i < ligne; i++)
    for (int j = 0; j < colonne; j++)
      t.m_values[t.idx(j, i)] = m_values[idx(i, j)];
  return t;
}

// ==========================================================================
// Accès aux éléments
// ==========================================================================

double &Matrix::operator()(int i, int j) {
  if (i < 0 || i >= ligne || j < 0 || j >= colonne)
    throw std::invalid_argument("Error: invalid index not in matrix");
  return m_values[idx(i, j)];
}

const double &Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= ligne || j < 0 || j >= colonne)
    throw std::invalid_argument("Error: invalid index not in matrix");
  return m_values[idx(i, j)];
}

std::vector<double> Matrix::getLigne(int index) const {
  if (index < 0 || index >= ligne)
    throw std::out_of_range("Error: invalid index <" + std::to_string(index) +
                            "> doesn't exist");
  auto debut = m_values.begin() + idx(index,0);
  return { debut, debut + colonne};
}

std::vector<double> Matrix::getCol(int index) const {
  if (index < 0 || index >= colonne)
    throw std::out_of_range("Error: invalid index <" + std::to_string(index) +
                            "> doesn't exist");
  std::vector<double> res;
  res.reserve(ligne);
  for (int i = 0; i < ligne; i++) {
    res.push_back(m_values[idx(i, index)]);
  }
  return res;
}

// ==========================================================================
// Opérateurs
// ==========================================================================

Matrix &Matrix::operator+=(const Matrix &other) {
  if (ligne != other.ligne || colonne != other.colonne)
    throw std::invalid_argument(
        "Error: Cannot add up 2 matrices of different size");
  for (std::size_t k = 0; k < m_values.size(); k++)
    m_values[k] += other.m_values[k];
  return *this;
}

Matrix &Matrix::operator-=(const Matrix &other) {
  if (ligne != other.ligne || colonne != other.colonne)
    throw std::invalid_argument(
        "Error: Cannot substract 2 matrices of different size");
  for (std::size_t k = 0; k < m_values.size(); k++)
    m_values[k] -= other.m_values[k];
  return *this;
}

Matrix operator+(Matrix lhs, const Matrix &rhs) { return lhs += rhs; }
Matrix operator-(Matrix lhs, const Matrix &rhs) { return lhs -= rhs; }

// carnet: debutProduit
Matrix Matrix::operator*(const Matrix &other) const {
  if (colonne != other.ligne)
    throw std::invalid_argument(
        "Error: cannot multiply, columns of A != lines of B");

  Matrix res(ligne, other.colonne); // rempli de 0
  const int p = other.colonne;
  for (int i = 0; i < ligne; i++) {
    double *ligne_res = &res.m_values[res.idx(i, 0)]; // ligne i du résultat
    for (int k = 0; k < colonne; k++) {
      const double a_ik = m_values[idx(i, k)];
      const double *ligne_b = &other.m_values[other.idx(k, 0)]; // ligne k de B
      for (int j = 0; j < p; j++)
        ligne_res[j] += a_ik * ligne_b[j];
    }
  }
  return res;
}
// carnet: finProduit

Matrix operator*(const Matrix &m, const std::vector<double> &vect) {
  if (m.colonne != static_cast<int>(vect.size()))
    throw std::invalid_argument(
        "Error: cannot multiply, vector size != columns of the matrix");

  Matrix res(m.ligne, 1);
  for (int i = 0; i < m.ligne; i++) {
    double somme = 0.0;
    for (int k = 0; k < m.colonne; k++)
      somme += m.m_values[m.idx(i, k)] * vect[k];
    res.m_values[i] = somme;
  }
  return res;
}

Matrix Matrix::operator+() const { return *this; }

Matrix Matrix::operator-() const {
  Matrix m(ligne, colonne);
  for (std::size_t k = 0; k < m_values.size(); k++)
    m.m_values[k] = -m_values[k];
  return m;
}

bool Matrix::operator==(const Matrix &m) const {
  if (ligne != m.ligne || colonne != m.colonne)
    throw std::invalid_argument(
        "Error: Cannot compare 2 matrices of different size");
  for (std::size_t k = 0; k < m_values.size(); k++)
    if (std::fabs(m_values[k] - m.m_values[k]) >= PRECISION)
      return false;
  return true;
}

// ==========================================================================
// Affichage
// ==========================================================================

std::ostream &operator<<(std::ostream &out, const Matrix &m) {
  const int precision = 3;

  // 1) Mesure : largeur du nombre le plus long, une fois formaté
  std::size_t largeur = 0;
  for (double x : m.m_values) {
    std::ostringstream t;
    t << std::fixed << std::setprecision(precision) << x;
    largeur = std::max(largeur, t.str().size());
  }

  // 2) Affichage, puis restauration de l'état du flux
  const auto precision_avant = out.precision();
  const auto flags_avant = out.flags();

  out << std::fixed << std::setprecision(precision) << "[[";
  for (int i = 0; i < m.ligne; i++) {
    if (i != 0)
      out << "]\n [";
    for (int j = 0; j < m.colonne; j++)
      out << ' ' << std::setw(static_cast<int>(largeur))
          << m.m_values[m.idx(i, j)];
  }
  out << "]]";

  out.precision(precision_avant);
  out.flags(flags_avant);
  return out;
}

} // namespace mnp
