#pragma once

#include <ostream>
#include <stdexcept>
#include <vector>

namespace mnp {

class Vector {
public:
  explicit Vector(std::vector<double> vec) : m_eles(vec) {}

  double norme() const;
  Vector normalisation() const;

  /// Operators
  Vector &operator+=(const Vector &x) {
    if (this->size() != x.size())
      throw std::invalid_argument(
          "Les deux vecteurs ne sont pas de la même taille");
    for (int i = 0; i < static_cast<int>(this->size()); i++) {
      m_eles[i] += x.m_eles.at(i);
    }
    return *this;
  }
  friend Vector operator+(Vector lhs, const Vector &rhs) {
    lhs += rhs;
    return lhs;
  }
  Vector operator-(const Vector &a) const;
  Vector operator*(const Vector &a) const;
  /// unary Operators
  Vector operator+() const;
  Vector operator-() const;
  
  friend Vector operator/(Vector x, const double &a) {
    for(int i = 0; i < static_cast<int>(x.size()) ; i++)
      x.m_eles[i] /=a;
    return x;
  }
  friend Vector operator/(const double &a, Vector x) { return x / a; }

  friend Vector operator*(Vector x, const double &a) {
    for (int i = 0; i < static_cast<int>(x.size()); i++)
      x.m_eles[i] *= a;
    return x;
  }
  friend Vector operator*(const double &a, Vector x) { return x * a; }

  bool operator==(const Vector &a) const = default;
  friend std::ostream &operator<<(std::ostream &os, const Vector &v) {
    os << "{";
    for (double eles : v.m_eles)
      os << eles << ",";

    os << "}";
    return os;
  }
  double &operator[](int index) {
    if (index < 0 || index > static_cast<int>(m_eles.size()))
      throw std::invalid_argument("L'index est introuvable dans le vecteur");
    return m_eles[index];
  }
  double operator[](int index) const {
    if (index < 0 || index > static_cast<int>(m_eles.size()))
      throw std::invalid_argument("L'index est introuvable dans le vecteur");
    return m_eles[index];
  }
  size_t size(void) const;

private:
  std::vector<double> m_eles;
};

} // namespace mnp
