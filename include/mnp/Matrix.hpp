#pragma once

#include <cstddef>
#define PRECISION 1e-9
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ios>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>
namespace mnp {
class Matrix {
public:
  // ========================================================================
  // Constructors
  // ========================================================================
  /**
   * @brief Create a matrix of size(i,j)
   *
   * @param i int : number of lines
   * @param j int : number of cols
   */
  explicit Matrix(const int i, const int j)
      : ligne(i), colonne(j), m_values(i * j, 0.0) {}

  /**
   * @brief Create a matrix from a vector of vector of double
   *
   * @param vect vector<vector<double>> values of the matrix
   */
  explicit Matrix(const std::vector<std::vector<double>> &vect)
      : ligne(vect.size()), colonne(vect.empty() ? 0 : vect[0].size()),
        m_values(colonne * ligne, 0.0) {
    if (colonne == 0 || ligne == 0) {
      throw std::invalid_argument("Vecteur de la mauvaise taille");
    }
    for (int i = 0; i < static_cast<int>(vect.size()); i++) {
      if (static_cast<int>(vect[i].size()) != colonne)
        throw std::invalid_argument("Vecteur de la mauvaise taille");
      for (int j = 0; j < static_cast<int>(vect[i].size()); j++) {
        m_values[(i * colonne) + j] = vect[i][j];
      }
    }
  }

  /**
   * @brief Generates an indentity matrix of size n
   *
   * @param n int size of matrix
   * @return Matrix square of size n
   */
  static Matrix identity(const int n) {
    Matrix m(n, n);
    for (int i = 0; i < (n * n); i++) {
      if ((i / n) == (i % n)) {
        m.m_values[i] = 1.0;
      }
    }
    return m;
  }

  /**
   * @brief Check if the matrix is square (size of : n * n )
   *
   * @return bool
   */
  bool isSquare() const { return ligne == colonne; }
  /**
   * @brief check if the matrix is filled by zeros
   *
   * @return bool
   */
  bool isZero() const {
    return std::all_of(m_values.begin(), m_values.end(),
                       [](double x) { return !(std::fabs(x) >= PRECISION); });
  }

  /**
   * @brief Check if the matrix is an identity matrix
   *
   * @return bool
   */
  bool isIndentity() const {
    // if not square it can't be an identity
    if (!isSquare())
      return false;

    for (int index = 0; index < static_cast<int>(m_values.size()); index++) {
      if (((index / colonne) == (index % colonne))) {
        if (std::fabs(1.0 - this->m_values[index]) >= PRECISION)
          return false;
      } else if (std::fabs(this->m_values[index]) >= PRECISION)
        return false;
    }
    return true;
  }

  /**
   * @brief transpose the Matrix from (i,j) -> (j,i)
   *
   * @return M^t
   */
  Matrix transpose() const {
    Matrix m(this->colonne, this->ligne);
    for (int i = 0; i < static_cast<int>(this->colonne * this->ligne); i++) {
      m.m_values[(i % this->colonne) * this->ligne + (i / this->colonne)] =
          this->m_values[i];
    }
    return m;
  }

  // ========================================================================
  // Operators
  // ========================================================================

  /**
   * @brief Get the value of (i,j) in the matrix
   *
   * @param i const int line of the matrix
   * @param j const int col of the matrix
   * @return value of i,j
   */
  double &operator()(const int i, const int j) {
    if (!(i < ligne && i >= 0) || !(j < colonne && j >= 0))
      throw std::invalid_argument("Error: invalid index not in matrix");
    return m_values.at(i * colonne + j);
  }

  /**
   * @brief Get the value of (i,j) in the matrix
   *
   * @param i const int line of the matrix
   * @param j const int col of the matrix
   * @return value of i,j
   */
  const double &operator()(const int i, const int j) const {
    if (!(i < ligne && i >= 0) || !(j < colonne && j >= 0))
      throw std::invalid_argument("Error: invalid index not in matrix");
    return m_values.at(i * colonne + j);
  }

  // binary operators
  Matrix &operator+=(const Matrix &other) {
    if (this->ligne != other.ligne || this->colonne != other.colonne)
      throw std::invalid_argument(
          "Error: Cannot add up 2 matrices of different size");
    for (int index = 0; index < static_cast<int>(m_values.size()); index++) {
      this->m_values[index] += other.m_values.at(index);
    }
    return *this;
  }

  Matrix &operator-=(const Matrix &other) {
    if (this->ligne != other.ligne || this->colonne != other.colonne)
      throw std::invalid_argument(
          "Error: Cannot substract 2 matrices of different size");
    for (int index = 0; index < static_cast<int>(m_values.size()); index++) {
      this->m_values[index] -= other.m_values.at(index);
    }
    return *this;
  }

  friend Matrix operator+(Matrix lhs, const Matrix &rhs) { return lhs += rhs; }
  friend Matrix operator-(Matrix lhs, const Matrix &rhs) { return lhs -= rhs; }

  /**
   * @brief Produit matriciel, ordre de boucles i, j, k (ancienne version)
   *
   * Calcule chaque case du résultat d'un coup : dans la boucle intérieure,
   * on lit B colonne par colonne (sauts de other.colonne en mémoire).
   * Gardée pour comparer les temps avec operator*.
   */
  Matrix produit_ijk(const Matrix &other) const {
    if (this->colonne != other.ligne)
      throw std::invalid_argument(
          "Error: cannot multiply, columns of A != lines of B");

    Matrix res(this->ligne, other.colonne);
    for (int index = 0; index < (res.ligne * res.colonne); index++) {
      double res_sum = 0;
      for (int i = 0; i < this->colonne; i++) {
        res_sum += (this->m_values[(i + this->colonne *
                                            (index / res.colonne))] // Ligne
                    * other.m_values[(i * other.colonne) +
                                     (index % res.colonne)] // Colonne
        );
      }
      res.m_values[index] = res_sum;
    }
    return res;
  }

  /**
   * @brief Produit matriciel, ordre de boucles i, k, j
   *
   * Pour chaque ligne i de A et chaque k, A(i,k) est lu une seule fois puis
   * multiplié par toute la ligne k de B, ajoutée à la ligne i du résultat.
   * Dans la boucle intérieure, B et le résultat sont parcourus case par case,
   * dans l'ordre où ils sont rangés en mémoire.
   */
  Matrix operator*(const Matrix &other) const {
    if (this->colonne != other.ligne)
      throw std::invalid_argument(
          "Error: cannot multiply, columns of A != lines of B");

    Matrix res(this->ligne, other.colonne); // rempli de 0
    for (int i = 0; i < this->ligne; i++) {
      for (int k = 0; k < this->colonne; k++) {
        const double a_ik = this->m_values[i * this->colonne + k];
        for (int j = 0; j < other.colonne; j++) {
          res.m_values[i * res.colonne + j] +=
              a_ik * other.m_values[k * other.colonne + j];
        }
      }
    }
    return res;
  }

  friend Matrix operator*(const Matrix &m, const std::vector<double> &vect) {
    if (m.colonne != static_cast<int>(vect.size()))
      throw std::invalid_argument(
          "Error: Cannot multiply a matrix by a vector if the size of the "
          "vector is not equal to the matrix j");

    Matrix res(m.ligne, 1);
    for (int i = 0; i < (res.ligne); i++) {

      double res_sum = 0;
      for (int index = 0; index < m.colonne; index++) {
        res_sum += (m.m_values[index + (m.colonne * (i / res.colonne))] // Ligne
                    * vect.at(index % static_cast<int>(vect.size())) // Vecteur
        );
      }
      res.m_values[i] = res_sum;
    }
    return res;
  }

  // unary Operators
  Matrix operator+() { return *this; }

  Matrix operator-() {
    Matrix m(ligne, colonne);
    for (int index = 0; index < static_cast<int>(m_values.size()); index++) {
      m.m_values[index] = this->m_values[index] * -1.0f;
    }
    return m;
  }

  bool operator==(const Matrix &m) const {
    if (this->ligne != m.ligne || this->colonne != m.colonne)
      throw std::invalid_argument(
          "Error: Cannot compare 2 matrices of different size");
    for (int i = 0; i < static_cast<int>(m_values.size()); i++)
      if (std::fabs(m_values.at(i) - m.m_values.at(i)) >= PRECISION)
        return false;
    return true;
  }

  // Affichage matrices
  friend std::ostream &operator<<(std::ostream &out, const Matrix &m) {
    std::size_t largeur = 0;
    const int precision = 3;
    for (double x : m.m_values) {
      std::ostringstream t;
      t << std::fixed << std::setprecision(precision) << x;
      largeur = std::max(largeur, t.str().size());
    }
    const auto default_precision(out.precision());
    const auto default_flags(out.flags());

    out << std::fixed << std::setprecision(precision) << "[[";
    for (int i = 0; i < m.ligne * m.colonne; i++) {
      if (i % m.colonne == 0 && i != 0) {
        out << "]\n [";
      }
      out << ' ' << std::setw(static_cast<int>(largeur)) << m.m_values[i];
    }
    out << "]]" << std::setprecision(default_precision);
    out.flags(default_flags);
    return out;
  }

  int nbLignes() const { return ligne; }
  int nbColonne() const { return colonne; }

private:
  const int ligne;
  const int colonne;
  std::vector<double> m_values;
};

} // namespace mnp
