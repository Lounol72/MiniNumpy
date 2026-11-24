#include "mnp/Vector.hpp"
#include "mnp/Precision.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace mnp {

// ==========================================================================
// Vérifications
// ==========================================================================

void Vector::verifierMemeTaille(const Vector &x) const {
  if (size() != x.size())
    throw std::invalid_argument("Les deux vecteurs ne sont pas de la même "
                                "taille (" +
                                std::to_string(size()) + " et " +
                                std::to_string(x.size()) + ")");
}

void Vector::verifierIndex(int index) const {
  if (index < 0 || static_cast<std::size_t>(index) >= size())
    throw std::out_of_range("L'index " + std::to_string(index) +
                            " n'est pas dans le vecteur (taille " +
                            std::to_string(size()) + ")");
}
// ==========================================================================
// Constructeurs
// ==========================================================================

Vector Vector::zeros(int n) {
  if (n < 0)
    throw std::invalid_argument(
        "Un vecteur doit au minimum avoir une taille 0");
  Vector v(std::vector<double>(n, 0.0));
  return v;
}
Vector Vector::ones(int n) {
  if (n < 0)
    throw std::invalid_argument(
        "Un vecteur doit au minimum avoir une taille 0");
  Vector v(std::vector<double>(n, 1.0));
  return v;
}

// ==========================================================================
// Propriétés
// ==========================================================================

double Vector::norme() const { return std::sqrt(scalaire_analytique(*this)); }

// carnet: debutNormalisation
Vector Vector::normalisation() const {
  const double n = norme();
  if (n == 0.0)
    throw std::invalid_argument("Normalisation impossible : vecteur nul");
  return *this / n;
}
// carnet: finNormalisation

double Vector::scalaire_analytique(const Vector &x) const {
  verifierMemeTaille(x);
  double res = 0.0;
  for (std::size_t i = 0; i < size(); i++)
    res += m_eles[i] * x.m_eles[i];
  return res;
}

double Vector::angle(const Vector &x) const {
  verifierMemeTaille(x);
  const double normes = norme() * x.norme();
  if (normes == 0.0)
    throw std::invalid_argument("Angle non défini avec un vecteur nul");
  const double cos_theta =
      std::clamp(scalaire_analytique(x) / normes, -1.0, 1.0);
  return std::acos(cos_theta);
}

double Vector::scalaire_geometrique(const Vector &x) const {
  verifierMemeTaille(x);
  const double normes = norme() * x.norme();
  if (normes == 0.0)
    return 0.0;
  return normes * std::cos(angle(x));
}

// ==========================================================================
// Accès aux éléments
// ==========================================================================

double &Vector::operator[](int index) {
  verifierIndex(index);
  return m_eles[index];
}

double Vector::operator[](int index) const {
  verifierIndex(index);
  return m_eles[index];
}

// ==========================================================================
// Opérateurs composés
// ==========================================================================

Vector &Vector::operator+=(const Vector &x) {
  verifierMemeTaille(x);
  for (std::size_t i = 0; i < size(); i++)
    m_eles[i] += x.m_eles[i];
  return *this;
}

Vector &Vector::operator-=(const Vector &x) {
  verifierMemeTaille(x);
  for (std::size_t i = 0; i < size(); i++)
    m_eles[i] -= x.m_eles[i];
  return *this;
}

Vector &Vector::operator*=(const Vector &x) {
  verifierMemeTaille(x);
  for (std::size_t i = 0; i < size(); i++)
    m_eles[i] *= x.m_eles[i];
  return *this;
}

Vector &Vector::operator*=(double a) {
  for (double &e : m_eles)
    e *= a;
  return *this;
}

Vector &Vector::operator/=(double a) {
  for (double &e : m_eles)
    e /= a;
  return *this;
}

// ==========================================================================
// Opérateurs binaires et unaires
// ==========================================================================

Vector operator/(double a, Vector x) {
  for (double &e : x.m_eles)
    e = a / e;
  return x;
}

Vector Vector::operator-() const { return *this * -1.0; }

bool Vector::operator==(const Vector &x) const {
  if (size() != x.size())
    return false;
  for (std::size_t i = 0; i < size(); i++)
    if (!approx_equal(m_eles[i], x.m_eles[i]))
      return false;
  return true;
}

std::ostream &operator<<(std::ostream &os, const Vector &v) {
  os << '{';
  for (std::size_t i = 0; i < v.size(); i++) {
    if (i != 0)
      os << ", ";
    os << v.m_eles[i];
  }
  return os << '}';
}

} // namespace mnp
