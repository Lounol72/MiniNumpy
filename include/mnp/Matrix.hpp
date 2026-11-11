#pragma once

#define PRECISION 1e-9

#include <algorithm>
#include <cmath>
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
  explicit Matrix(std::vector<std::vector<double>> &vect)
      : ligne(vect.size()), colonne(vect[0].size()),
        m_values(colonne * ligne, 0.0) {
    for (int i = 0; i < static_cast<int>(vect.size()); i++)
      for (int j = 0; j < static_cast<int>(vect[i].size()); j++) {
        m_values[(i * colonne) + j] = vect[i][j];
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

  Matrix operator*(const Matrix &other) const {
    if (this->colonne != other.ligne)
      throw std::invalid_argument(
          "Error: Cannot substract 2 matrices of different size");

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
                    * vect.at(index % static_cast<int>(vect.size()))  // Vecteur
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
  int nbLignes() const { return ligne; }
  int nbColonne() const { return colonne; }

private:
  const int ligne;
  const int colonne;
  std::vector<double> m_values;
};

} // namespace mnp
