#pragma once

#include <algorithm>
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
      for (int j = 0; j < static_cast<int>(vect.size()); j++) {
        m_values[(i * ligne) + j] = vect[i][j];
      }
  }

  bool isSquare() { return ligne == colonne; }
  bool isZero() {
    return std::all_of(m_values.begin(), m_values.end(),
                       [](int x) { return x == 0; });
  }

  bool isIndentity() {
    if (!isSquare())
      return false;
    if (m_values[0] != 1)
      return false;

    for (int index = 1; index < static_cast<int>(m_values.size()); index++) {
      if (((index / ligne) == (index % ligne)) && m_values[index] != 1)
        return false;
      else if (m_values[index] != 0)
        return false;
    }
    return true;
  }
  
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
    for (int index = 0; index < (this->ligne * other.colonne); index++) {
      double res_sum = 0;
      for (int i = 0; i < this->ligne; i++) {
        res_sum +=
            (this->m_values[(i * res.ligne) + (index % res.ligne)]   // Colonne
             * other.m_values[i + (res.ligne * (index / res.ligne))] // Ligne
            );
      }
      res.m_values[index] = res_sum;
    }
    return res;
  }

  // unary Operators
  Matrix &operator+() {
    for (int index = 0; index < static_cast<int>(m_values.size()); index++) {
      this->m_values[index] = this->m_values[index];
    }
    return *this;
  }

  Matrix &operator-() {
    for (int index = 0; index < static_cast<int>(m_values.size()); index++) {
      this->m_values[index] = this->m_values[index] * -1;
    }
    return *this;
  }

  bool operator==(const Matrix &m) const {
    if (this->ligne != m.ligne || this->colonne != m.colonne)
      throw std::invalid_argument(
          "Error: Cannot compare 2 matrices of different size");
    for (int i = 0; i < static_cast<int>(m_values.size()); i++)
      if (m_values.at(i) != m.m_values.at(i))
        return false;
    return true;
  }
  bool nbLignes() const { return ligne; }
  bool nbColonne() const { return colonne; }

private:
  const int ligne;
  const int colonne;
  std::vector<double> m_values;
};

} // namespace mnp
