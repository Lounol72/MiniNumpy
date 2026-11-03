#include <stdexcept>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "mnp/version.hpp"
#include "mnp/Vector.hpp"

TEST_CASE("vérification de version"){
  CHECK_MESSAGE(mnp::version()=="0.0.1", "La version doit suivre le tag git");
}

TEST_CASE("Vérifications des vecteurs opérations + et +="){
  mnp::Vector u({3,0,4});
  mnp::Vector v({1,2,2});
  mnp::Vector w({1,2});
  CHECK_MESSAGE((u + v) == mnp::Vector({4,2,6}),"L'addition vectorielle n'est pas correcte");
  CHECK_THROWS_AS(u + w, std::invalid_argument);
  //CHECK_MESSAGE((u * v) == mnp::Vector({3,0,8}),"La multiplication vectorielle n'est pas correcte");
  //CHECK_MESSAGE((2.0 * u) == mnp::Vector({6,0,8}),"La multiplication vectorielle n'est pas correcte");
}
