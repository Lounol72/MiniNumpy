#include "mnp/Vector.hpp"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace mnp {

double Vector::norme() const {
  return std::sqrt(this->scalaire_analytique(*this));
}
Vector Vector::normalisation() const {
  double norme = this->norme();
  if (norme == 0.0)
    throw std::invalid_argument("Normalisation impossible norme =" +
                                std::to_string(norme));
  std::vector<double> res;

  for (int i = 0; i < static_cast<int>(this->size()); i++) {
    res.push_back(this->m_eles.at(i) / norme);
  }
  return Vector(res);
}

Vector Vector::operator-(const Vector &a) const {
  if (this->size() != a.size())
    throw std::invalid_argument(
        "Les deux vecteurs ne sont pas de la même taille");
  std::vector<double> res;

  for (int i = 0; i < static_cast<int>(this->size()); i++) {
    res.push_back(this->m_eles.at(i) - a.m_eles.at(i));
  }
  return Vector(res);
}
Vector Vector::operator*(const Vector &a) const {
  if (this->size() != a.size())
    throw std::invalid_argument(
        "Les deux vecteurs ne sont pas de la même taille");
  std::vector<double> res;

  for (int i = 0; i < static_cast<int>(this->size()); i++) {
    res.push_back(this->m_eles.at(i) * a.m_eles.at(i));
  }
  return Vector(res);
}
/// unary Operators
Vector Vector::operator+() const {
  std::vector<double> res;

  for (int i = 0; i < static_cast<int>(this->size()); i++) {
    res.push_back(this->m_eles.at(i));
  }
  return Vector(res);
}
Vector Vector::operator-() const {
  std::vector<double> res;

  for (int i = 0; i < static_cast<int>(this->size()); i++) {
    res.push_back(this->m_eles.at(i) * -1);
  }
  return Vector(res);
}

size_t Vector::size(void) const { return m_eles.size(); }
} // namespace mnp
