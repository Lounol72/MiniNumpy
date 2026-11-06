#include <cmath>
#include <numbers>
#include <stdexcept>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "mnp/Vector.hpp"
#include "mnp/version.hpp"

TEST_CASE("vérification de version") {
  CHECK_MESSAGE(mnp::version() == "0.0.1", "La version doit suivre le tag git");
}

TEST_CASE("Vérifications des opérations vectorielles :  +, +=, -= et -") {
  mnp::Vector u({3, 0, 4});
  mnp::Vector v({1, 2, 2});
  mnp::Vector w({1, 2});
  CHECK_MESSAGE((u + v) == mnp::Vector({4, 2, 6}),
                "L'addition vectorielle n'est pas correcte");
  CHECK_MESSAGE((u + v) == (v + u),
                "L'addition vectorielle est sensée être commutative");
  CHECK_THROWS_AS(u + w, std::invalid_argument);
  mnp::Vector somme((u + v));
  u += v;
  CHECK_MESSAGE(somme == u, "L'addition vectorielle n'est pas correcte");
  CHECK_MESSAGE((u -= v) == mnp::Vector({3, 0, 4}),
                "La soustraction vectorielle n'est pas correcte");
  CHECK_MESSAGE((u += v) == mnp::Vector({4, 2, 6}),
                "L'addition vectorielle n'est pas correcte");
  CHECK_MESSAGE((u -= v) == mnp::Vector({3, 0, 4}),
                "La soustraction vectorielle n'est pas correcte");
  CHECK_MESSAGE((u - v) == mnp::Vector({2, -2, 2}),
                "La soustraction vectorielle n'est pas correcte");
}

TEST_CASE("Vérifications des opérations vectorielles : produit scalaire, "
          "produit par un scalaire, norme, normalisation") {
  mnp::Vector u({3, 0, 4});
  mnp::Vector v({1, 2, 2});
  mnp::Vector w({1, 2});
  CHECK_MESSAGE((u * v) == mnp::Vector({3, 0, 8}),
                "Le produit scalaire n'est pas correct");
  CHECK_MESSAGE((u * v) == (v * u), "Le produit scalaire n'est pas correct");
  CHECK_MESSAGE((2.0 * u) == mnp::Vector({6, 0, 8}),
                "La produit par scalaire n'est pas correct");
  CHECK_MESSAGE((u * 2.0) == (2.0 * u),
                "La produit par scalaire n'est pas correct");
  CHECK_MESSAGE((u.norme()) == 5.0f, "La norme n'est pas correct");
  CHECK_MESSAGE((u.norme()) == std::sqrt(u.scalaire_analytique(u)),
                "La norme n'est pas correct");
  CHECK_MESSAGE((u.normalisation()) == mnp::Vector({0.6, 0, 0.8}),
                "La normalisation n'est pas correcte");
  CHECK_MESSAGE((u.normalisation().norme()) == 1.0f,
                "La normalisation n'est pas correct");
  double normes_mult = u.norme() * v.norme();
  double angle = u.scalaire_analytique(v) / (normes_mult);
  CHECK_MESSAGE((normes_mult * std::acos(angle)) ==
                    doctest::Approx(11.21383518f).epsilon(1e-5),
                "L'angle n'est pas correct");
  CHECK_MESSAGE(((std::acos(angle)) >= 0.0f),
                "L'angle n'est pas correct");
   CHECK_MESSAGE((( std::acos(angle)) <= std::numbers::pi),
                "L'angle n'est pas correct");
}
