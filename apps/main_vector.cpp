// Étape 1 : angle entre deux vecteurs avec mnp::Vector.
// À comparer avec le calcul à la main de la section 1 du guide.
#include <cmath>
#include <iostream>
#include <numbers>

#include "mnp/Vector.hpp"

int main() {
  mnp::Vector u({1, 2, 2});
  mnp::Vector v({3, 0, 4});

  double produit = u.scalaire_analytique(v);
  double cos_theta = produit / (u.norme() * v.norme());
  double theta = std::acos(cos_theta);

  std::cout << "u = " << u << "\n";
  std::cout << "v = " << v << "\n";
  std::cout << "u.v = " << produit << "\n";
  std::cout << "|u| = " << u.norme() << ", |v| = " << v.norme() << "\n";
  std::cout << "cos(theta) = " << cos_theta << "\n";
  std::cout << "theta = " << theta << " rad = "
            << theta * 180.0 / std::numbers::pi << " deg\n";
  return 0;
}
