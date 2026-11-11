#include <cmath>
#include <numbers>
#include <sstream>
#include <stdexcept>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "mnp/Matrix.hpp"
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
  CHECK_MESSAGE(((std::acos(angle)) >= 0.0f), "L'angle n'est pas correct");
  CHECK_MESSAGE(((std::acos(angle)) <= std::numbers::pi),
                "L'angle n'est pas correct");
}

TEST_CASE("Vecteur : taille et accès [i]") {
  mnp::Vector u({3, 0, 4});
  CHECK(u.size() == 3);
  CHECK(u[0] == 3);
  CHECK(u[2] == 4);
  u[1] = 7;
  CHECK_MESSAGE(u == mnp::Vector({3, 7, 4}), "u[1] = 7 doit modifier u");
  CHECK_THROWS_AS(u[-1], std::invalid_argument);
  CHECK_THROWS_AS(u[3], std::invalid_argument);

  const mnp::Vector cu({1, 2});
  CHECK(cu[1] == 2);
  CHECK_THROWS_AS(cu[2], std::invalid_argument);
}

TEST_CASE("Vecteur : tailles incompatibles") {
  mnp::Vector u({3, 0, 4});
  mnp::Vector w({1, 2});
  CHECK_THROWS_AS(u - w, std::invalid_argument);
  CHECK_THROWS_AS(u * w, std::invalid_argument);
  CHECK_THROWS_AS(u += w, std::invalid_argument);
  CHECK_THROWS_AS(u -= w, std::invalid_argument);
  CHECK_THROWS_AS((void)(u == w), std::invalid_argument);
  CHECK_THROWS_AS(u.scalaire_analytique(w), std::invalid_argument);
  CHECK_THROWS_AS(u.scalaire_geometrique(w), std::invalid_argument);
}

TEST_CASE("Vecteur : opérateurs unaires") {
  mnp::Vector u({3, -1, 0});
  CHECK(+u == mnp::Vector({3, -1, 0}));
  CHECK(-u == mnp::Vector({-3, 1, 0}));
  CHECK_MESSAGE(-(-u) == u, "-(-u) doit valoir u");
  CHECK_MESSAGE(u == mnp::Vector({3, -1, 0}), "-u ne doit pas modifier u");
  CHECK(u + (-u) == mnp::Vector({0, 0, 0}));
}

TEST_CASE("Vecteur : division par un scalaire") {
  mnp::Vector u({3, 0, 4});
  CHECK(u / 2.0 == mnp::Vector({1.5, 0, 2}));
  CHECK_MESSAGE((u * 2.0) / 2.0 == u, "(u * 2) / 2 doit valoir u");
  mnp::Vector v({1, 2, 4});
  CHECK_MESSAGE(
      4.0 / v == mnp::Vector({4, 2, 1}),
      "a / v doit calculer a / v_i pour chaque élément (comme NumPy)");
}

TEST_CASE("Vecteur : produit scalaire analytique") {
  mnp::Vector u({3, 0, 4});
  mnp::Vector v({1, 2, 2});
  CHECK(u.scalaire_analytique(v) == doctest::Approx(11.0));
  CHECK_MESSAGE(u.scalaire_analytique(v) == v.scalaire_analytique(u),
                "Le produit scalaire doit être symétrique");
  CHECK(u.scalaire_analytique(u) == doctest::Approx(25.0));
  mnp::Vector e1({1, 0});
  mnp::Vector e2({0, 1});
  CHECK_MESSAGE(e1.scalaire_analytique(e2) == doctest::Approx(0.0),
                "Deux vecteurs orthogonaux ont un produit scalaire nul");
}

TEST_CASE("Vecteur : produit scalaire géométrique") {
  mnp::Vector u({3, 0, 4});
  mnp::Vector v({1, 2, 2});
  CHECK_MESSAGE(u.scalaire_geometrique(v) ==
                    doctest::Approx(u.scalaire_analytique(v)),
                "||u||·||v||·cos(θ) doit donner le même résultat que la "
                "formule analytique");
  mnp::Vector e1({1, 0});
  mnp::Vector e2({0, 1});
  CHECK(e1.scalaire_geometrique(e2) == doctest::Approx(0.0));
}

TEST_CASE("Vecteur : norme et normalisation, cas limites") {
  mnp::Vector nul({0, 0, 0});
  CHECK(nul.norme() == 0.0);
  CHECK_THROWS_AS(nul.normalisation(), std::invalid_argument);
  mnp::Vector u({-3, 0, -4});
  CHECK_MESSAGE(u.norme() == doctest::Approx(5.0),
                "La norme ne dépend pas du signe");
  CHECK(u.normalisation() == mnp::Vector({-0.6, 0, -0.8}));
}

TEST_CASE("Vecteur : affichage <<") {
  std::ostringstream os;
  os << mnp::Vector({1, 2.5, -3});
  CHECK(os.str() == "{1,2.5,-3,}");
}

TEST_CASE("Vérifications des matrices") {
  mnp::Matrix a(3, 3);
  mnp::Matrix b(1, 3);
  mnp::Matrix c(3, 3);
  SUBCASE("Remplissage des matrices") {
    a(0, 0) = 1;
    c(0, 0) = 1;
    CHECK_MESSAGE(a == c, "Le remplissage de matrices ne fonctionnes pas");
    CHECK_THROWS_AS(a == b, std::invalid_argument);
  }
  SUBCASE("Produit Matriciel") {
    for (int i = 0; i < 3 * 3; i++) {
      a(i % a.nbLignes(), i % a.nbColonne()) = i;
      c(i % c.nbLignes(), i % c.nbColonne()) = i;
    }
  }
}

// ============================================================================
// Tests des matrices
// ============================================================================

namespace {
// Construit une matrice élément par élément via operator(), pour que les tests
// ne dépendent pas du constructeur vector<vector<double>>.
mnp::Matrix fabriquer(const std::vector<std::vector<double>> &valeurs) {
  const int lignes = static_cast<int>(valeurs.size());
  const int colonnes = static_cast<int>(valeurs[0].size());
  mnp::Matrix m(lignes, colonnes);
  for (int i = 0; i < lignes; i++)
    for (int j = 0; j < colonnes; j++)
      m(i, j) = valeurs[i][j];
  return m;
}
} // namespace

TEST_CASE("Matrice : constructeur (lignes, colonnes)") {
  mnp::Matrix m(2, 3);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++)
      CHECK_MESSAGE(m(i, j) == 0.0, "Une nouvelle matrice doit être nulle");
}

TEST_CASE("Matrice : dimensions") {
  mnp::Matrix m(2, 3);
  const int lignes = 2;
  const int colonnes = 3;
  CHECK_MESSAGE(m.nbLignes() == lignes, "nbLignes() doit renvoyer 2");
  CHECK_MESSAGE(m.nbColonne() == colonnes, "nbColonne() doit renvoyer 3");
}

TEST_CASE("Matrice : constructeur vector<vector<double>> (carrée)") {
  std::vector<std::vector<double>> v = {{1, 2}, {3, 4}};
  mnp::Matrix m(v);
  CHECK(m(0, 0) == 1);
  CHECK(m(0, 1) == 2);
  CHECK(m(1, 0) == 3);
  CHECK(m(1, 1) == 4);
}

TEST_CASE("Matrice : constructeur vector<vector<double>> (2x3)") {
  // NB : pas de test 3x2 ici, il provoque actuellement une lecture hors
  // limites (ASan arrête alors tout l'exécutable).
  std::vector<std::vector<double>> v = {{1, 2, 3}, {4, 5, 6}};
  mnp::Matrix m(v);
  CHECK_MESSAGE(
      m == fabriquer({{1, 2, 3}, {4, 5, 6}}),
      "Le constructeur doit recopier toutes les valeurs à leur place");
}

TEST_CASE("Matrice : accès (i, j) et bornes") {
  mnp::Matrix m(2, 3);
  m(1, 2) = 7.5;
  CHECK(m(1, 2) == 7.5);
  CHECK_MESSAGE(m(0, 0) == 0.0,
                "Écrire en (1,2) ne doit rien changer ailleurs");

  CHECK_THROWS_AS(m(-1, 0), std::invalid_argument);
  CHECK_THROWS_AS(m(0, -1), std::invalid_argument);
  CHECK_THROWS_AS(m(2, 0), std::invalid_argument);
  CHECK_THROWS_AS(m(0, 3), std::invalid_argument);

  const mnp::Matrix &cm = m;
  CHECK(cm(1, 2) == 7.5);
  CHECK_THROWS_AS(cm(2, 3), std::invalid_argument);
}

TEST_CASE("Matrice : égalité") {
  mnp::Matrix a = fabriquer({{1, 2}, {3, 4}});
  mnp::Matrix b = fabriquer({{1, 2}, {3, 4}});
  mnp::Matrix c = fabriquer({{1, 2}, {3, 5}});
  CHECK(a == b);
  CHECK_FALSE(a == c);
  b(0, 0) += 1e-12;
  CHECK_MESSAGE(a == b, "L'égalité doit tolérer un écart < PRECISION");
  CHECK_THROWS_AS(a == mnp::Matrix(2, 3), std::invalid_argument);
}

TEST_CASE("Matrice : identité, isSquare, isIndentity") {
  mnp::Matrix i3 = mnp::Matrix::identity(3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      CHECK(i3(i, j) == (i == j ? 1.0 : 0.0));

  CHECK(i3.isSquare());
  CHECK_FALSE(mnp::Matrix(2, 3).isSquare());

  CHECK(i3.isIndentity());
  CHECK_FALSE(mnp::Matrix(3, 3).isIndentity());
  CHECK_FALSE(mnp::Matrix(2, 3).isIndentity());
  mnp::Matrix presque = mnp::Matrix::identity(3);
  presque(0, 2) = 0.5;
  CHECK_FALSE(presque.isIndentity());
}

TEST_CASE("Matrice : isZero") {
  CHECK_MESSAGE(mnp::Matrix(2, 3).isZero(),
                "Une matrice remplie de 0 doit être nulle");
  mnp::Matrix m(2, 2);
  m(1, 1) = 1.0;
  CHECK_FALSE(m.isZero());
  m(1, 1) = 0.5;
  CHECK_FALSE_MESSAGE(m.isZero(),
                      "0.5 n'est pas nul (attention aux conversions)");
}

TEST_CASE("Matrice : transposée") {
  mnp::Matrix a = fabriquer({{1, 2, 3}, {4, 5, 6}});
  mnp::Matrix at = a.transpose();
  CHECK(at == fabriquer({{1, 4}, {2, 5}, {3, 6}}));
  CHECK_MESSAGE(at.transpose() == a, "(A^T)^T doit valoir A");
  mnp::Matrix i3 = mnp::Matrix::identity(3);
  CHECK(i3.transpose() == i3);
}

TEST_CASE("Matrice : addition et soustraction") {
  mnp::Matrix a = fabriquer({{1, 2, 3}, {4, 5, 6}});
  mnp::Matrix b = fabriquer({{6, 5, 4}, {3, 2, 1}});

  CHECK(a + b == fabriquer({{7, 7, 7}, {7, 7, 7}}));
  CHECK_MESSAGE(a + b == b + a, "L'addition doit être commutative");
  CHECK(a - b == fabriquer({{-5, -3, -1}, {1, 3, 5}}));
  CHECK(a - a == mnp::Matrix(2, 3));
  CHECK_MESSAGE(a == fabriquer({{1, 2, 3}, {4, 5, 6}}),
                "a + b et a - b ne doivent pas modifier a");

  mnp::Matrix c = a;
  c += b;
  CHECK(c == a + b);
  c -= b;
  CHECK(c == a);

  mnp::Matrix autre(3, 2);
  CHECK_THROWS_AS(a + autre, std::invalid_argument);
  CHECK_THROWS_AS(a - autre, std::invalid_argument);
  CHECK_THROWS_AS(c += autre, std::invalid_argument);
  CHECK_THROWS_AS(c -= autre, std::invalid_argument);
}

TEST_CASE("Matrice : opérateurs unaires") {
  mnp::Matrix a = fabriquer({{1, -2}, {3, 0}});
  mnp::Matrix plus = +a;
  CHECK(plus == fabriquer({{1, -2}, {3, 0}}));
  mnp::Matrix moins = -a;
  CHECK(moins == fabriquer({{-1, 2}, {-3, 0}}));
  CHECK_MESSAGE(a == fabriquer({{1, -2}, {3, 0}}),
                "-a ne doit pas modifier a (comme pour un double)");
}

TEST_CASE("Matrice : produit matriciel (carré)") {
  mnp::Matrix a = fabriquer({{1, 2}, {3, 4}});
  mnp::Matrix b = fabriquer({{5, 6}, {7, 8}});
  CHECK_MESSAGE(a * b == fabriquer({{19, 22}, {43, 50}}),
                "A*B : ligne i de A · colonne j de B");
  CHECK_MESSAGE(b * a == fabriquer({{23, 34}, {31, 46}}),
                "B*A : le produit n'est pas commutatif");
  mnp::Matrix i2 = mnp::Matrix::identity(2);
  CHECK(a * i2 == a);
  CHECK(i2 * a == a);
}

TEST_CASE("Matrice : produit matriciel (rectangulaire)") {
  // NB : pas de test (3x2)*(2x3) ici, il provoque actuellement une lecture
  // hors limites (ASan arrête alors tout l'exécutable).
  mnp::Matrix a = fabriquer({{1, 2, 3}, {4, 5, 6}});      // 2x3
  mnp::Matrix b = fabriquer({{7, 8}, {9, 10}, {11, 12}}); // 3x2
  CHECK_MESSAGE(a * b == fabriquer({{58, 64}, {139, 154}}),
                "(2x3)*(3x2) doit donner une 2x2");

  CHECK_MESSAGE(a * mnp::Matrix::identity(3) == a, "A * I3 doit valoir A");
  CHECK_MESSAGE(mnp::Matrix::identity(2) * a == a, "I2 * A doit valoir A");

  mnp::Matrix ligne = fabriquer({{1, 2, 3}});       // 1x3
  mnp::Matrix colonne = fabriquer({{4}, {5}, {6}}); // 3x1
  CHECK_MESSAGE(ligne * colonne == fabriquer({{32}}),
                "(1x3)*(3x1) = produit scalaire");

  CHECK_THROWS_AS(a * a, std::invalid_argument);
}

TEST_CASE("Matrice : produit matrice × vecteur (carré)") {
  mnp::Matrix a = fabriquer({{1, 2}, {3, 4}});
  std::vector<double> x = {5, 6};
  CHECK(a * x == fabriquer({{17}, {39}}));
  CHECK(mnp::Matrix::identity(2) * x == fabriquer({{5}, {6}}));
}

TEST_CASE("Matrice : produit matrice × vecteur (rectangulaire)") {
  // NB : pas de test avec plus de lignes que de colonnes ici, il provoque
  // actuellement une écriture hors limites (ASan arrête tout l'exécutable).
  mnp::Matrix a = fabriquer({{1, 2, 3}, {4, 5, 6}}); // 2x3
  std::vector<double> x = {1, 0, -1};
  CHECK_MESSAGE(a * x == fabriquer({{-2}, {-2}}),
                "(2x3) * vecteur de taille 3 doit donner une 2x1");

  std::vector<double> mauvais = {1, 2};
  CHECK_THROWS_AS(a * mauvais, std::invalid_argument);
}
