#include "mnp/Linalg.hpp"
#include <cmath>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <vector>
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
  std::vector<std::vector<double>> v = {{1, 2, 3}, {4, 5, 6}};
  mnp::Matrix m(v);
  CHECK_MESSAGE(
      m == fabriquer({{1, 2, 3}, {4, 5, 6}}),
      "Le constructeur doit recopier toutes les valeurs à leur place");
}

TEST_CASE("Matrice : constructeur vector<vector<double>> (3x2)") {
  mnp::Matrix m({{1, 2}, {3, 4}, {5, 6}});
  CHECK(m.nbLignes() == 3);
  CHECK(m.nbColonne() == 2);
  CHECK(m == fabriquer({{1, 2}, {3, 4}, {5, 6}}));
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
  mnp::Matrix a = fabriquer({{1, 2, 3}, {4, 5, 6}});      // 2x3
  mnp::Matrix b = fabriquer({{7, 8}, {9, 10}, {11, 12}}); // 3x2
  CHECK_MESSAGE(a * b == fabriquer({{58, 64}, {139, 154}}),
                "(2x3)*(3x2) doit donner une 2x2");
  CHECK_MESSAGE(b * a == fabriquer({{39, 54, 69}, {49, 68, 87}, {59, 82, 105}}),
                "(3x2)*(2x3) doit donner une 3x3");

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
  mnp::Matrix a = fabriquer({{1, 2, 3}, {4, 5, 6}}); // 2x3
  std::vector<double> x = {1, 0, -1};
  CHECK_MESSAGE(a * x == fabriquer({{-2}, {-2}}),
                "(2x3) * vecteur de taille 3 doit donner une 2x1");

  mnp::Matrix haute = fabriquer({{1, 2}, {3, 4}, {5, 6}}); // 3x2
  std::vector<double> y = {1, -1};
  CHECK_MESSAGE(haute * y == fabriquer({{-1}, {-1}, {-1}}),
                "(3x2) * vecteur de taille 2 doit donner une 3x1");

  std::vector<double> mauvais = {1, 2};
  CHECK_THROWS_AS(a * mauvais, std::invalid_argument);
}

TEST_CASE("Matrice : affichage <<") {
  SUBCASE("colonnes alignées sur le nombre le plus long") {
    // largeur = 7 ("333.000" et "-22.500"), plus l'espace devant chaque nombre
    std::ostringstream os;
    os << fabriquer({{1, -22.5}, {333, 0}});
    CHECK(os.str() == "[[   1.000 -22.500]\n"
                      " [ 333.000   0.000]]");
  }

  SUBCASE("matrice 2x3 de même largeur") {
    std::ostringstream os;
    os << fabriquer({{1, 2, 3}, {4, 5, 6}});
    CHECK(os.str() == "[[ 1.000 2.000 3.000]\n"
                      " [ 4.000 5.000 6.000]]");
  }

  SUBCASE("cas limites : 1x1 et vecteur colonne") {
    std::ostringstream os;
    os << fabriquer({{7}});
    CHECK(os.str() == "[[ 7.000]]");

    std::ostringstream col;
    col << fabriquer({{1}, {-10}});
    CHECK(col.str() == "[[   1.000]\n"
                       " [ -10.000]]");
  }

  SUBCASE("l'état du flux est restauré après l'affichage") {
    std::ostringstream os;
    const auto precision_avant = os.precision();
    const auto flags_avant = os.flags();
    os << fabriquer({{1, 2}, {3, 4}});
    CHECK_MESSAGE(os.precision() == precision_avant,
                  "operator<< doit restaurer la précision du flux");
    CHECK_MESSAGE(os.flags() == flags_avant,
                  "operator<< doit restaurer les flags (fixed...)");

    std::ostringstream pi;
    pi << fabriquer({{1}}) << '|' << 3.14159265358979;
    CHECK_MESSAGE(pi.str() == "[[ 1.000]]|3.14159",
                  "un double écrit après une matrice garde le format par "
                  "défaut");
  }
}

TEST_CASE("Matrice : produit, comparaison avec un calcul de référence") {
  // Référence volontairement naïve, écrite avec operator() : C(i,j) est la
  // somme des A(i,k) * B(k,j). Elle ne partage aucun code avec operator*.
  auto reference = [](const mnp::Matrix &a, const mnp::Matrix &b) {
    mnp::Matrix c(a.nbLignes(), b.nbColonne());
    for (int i = 0; i < a.nbLignes(); i++)
      for (int j = 0; j < b.nbColonne(); j++)
        for (int k = 0; k < a.nbColonne(); k++)
          c(i, j) += a(i, k) * b(k, j);
    return c;
  };

  SUBCASE("non carrées : (2x3) * (3x4)") {
    mnp::Matrix a = fabriquer({{1, -2, 3}, {0.5, 4, -1}});
    mnp::Matrix b = fabriquer({{2, 0, 1, -3}, {1, 1, 0, 2}, {-1, 5, 2, 0.25}});
    CHECK(a * b == fabriquer({{-3, 13, 7, -6.25}, {6, -1, -1.5, 6.25}}));
    CHECK(a * b == reference(a, b));
  }

  SUBCASE("rectangulaires plus grandes : (7x5) * (5x9)") {
    mnp::Matrix a(7, 5);
    mnp::Matrix b(5, 9);
    for (int i = 0; i < 7; i++)
      for (int j = 0; j < 5; j++)
        a(i, j) = (i * 3 + j * 7) % 11 - 5;
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 9; j++)
        b(i, j) = (i * 5 + j * 2) % 13 - 6;
    CHECK(a * b == reference(a, b));
  }

  SUBCASE("dimensions incompatibles") {
    mnp::Matrix a(2, 3);
    CHECK_THROWS_AS(a * a, std::invalid_argument);
  }
}

TEST_CASE("Matrice : constructeur vector<vector<double>>, entrées invalides") {
  SUBCASE("lignes de longueurs différentes") {
    CHECK_THROWS_AS(mnp::Matrix({{1, 2}, {3}}), std::invalid_argument);
    CHECK_THROWS_AS(mnp::Matrix({{1}, {2, 3}}), std::invalid_argument);
    CHECK_THROWS_AS(mnp::Matrix({{1, 2, 3}, {4, 5, 6}, {7, 8}}),
                    std::invalid_argument);
  }
  SUBCASE("matrice vide") {
    CHECK_THROWS_AS(mnp::Matrix(std::vector<std::vector<double>>{}),
                    std::invalid_argument);
    CHECK_THROWS_AS(mnp::Matrix({{}}), std::invalid_argument);
  }
}

TEST_CASE("Matrice : constructeur depuis une liste écrite directement") {
  mnp::Matrix m({{1, 2, 3}, {4, 5, 6}}); // possible grâce au const &
  CHECK(m.nbLignes() == 2);
  CHECK(m.nbColonne() == 3);
  CHECK(m == fabriquer({{1, 2, 3}, {4, 5, 6}}));

  mnp::Matrix colonne({{1}, {2}, {3}}); // 3x1 : l'ancien cas « NB »
  CHECK(colonne.nbLignes() == 3);
  CHECK(colonne.nbColonne() == 1);
  CHECK(colonne(2, 0) == 3);
}

// ============================================================================
// Calculs à la main de la note (étape 2, section 1)
// A = [[1, 2], [3, 4]], B = [[0, 1], [1, 0]], C = [[1, 2, 3], [4, 5, 6]]
// ============================================================================

TEST_CASE("Matrice : calculs à la main de la note") {
  const mnp::Matrix A({{1, 2}, {3, 4}});
  const mnp::Matrix B({{0, 1}, {1, 0}});
  const mnp::Matrix C({{1, 2, 3}, {4, 5, 6}});

  SUBCASE(
      "1. AB et BA : B à droite échange les colonnes, à gauche les lignes") {
    CHECK(A * B == mnp::Matrix({{2, 1}, {4, 3}}));
    CHECK(B * A == mnp::Matrix({{3, 4}, {1, 2}}));
    CHECK_FALSE_MESSAGE(A * B == B * A, "Le produit n'est pas commutatif");
  }

  SUBCASE("2. Transposées et tailles des produits") {
    CHECK(A.transpose() == mnp::Matrix({{1, 3}, {2, 4}}));
    CHECK(C.transpose() == mnp::Matrix({{1, 4}, {2, 5}, {3, 6}}));

    const mnp::Matrix Ct = C.transpose();
    const mnp::Matrix AC = A * C; // (2x2)(2x3) -> 2x3
    CHECK(AC.nbLignes() == 2);
    CHECK(AC.nbColonne() == 3);

    const mnp::Matrix CtA = Ct * A; // (3x2)(2x2) -> 3x2
    CHECK(CtA == mnp::Matrix({{13, 18}, {17, 24}, {21, 30}}));

    const mnp::Matrix CCt = C * Ct; // (2x3)(3x2) -> 2x2
    CHECK(CCt == mnp::Matrix({{14, 32}, {32, 77}}));

    const mnp::Matrix CtC = Ct * C; // (3x2)(2x3) -> 3x3
    CHECK(CtC == mnp::Matrix({{17, 22, 27}, {22, 29, 36}, {27, 36, 45}}));

    // CA : (2x3)(2x2), 3 colonnes != 2 lignes -> n'existe pas
    CHECK_THROWS_AS(C * A, std::invalid_argument);
  }

  SUBCASE("3. AC, puis (AC)^T = C^T A^T") {
    const mnp::Matrix AC = A * C;
    CHECK(AC == mnp::Matrix({{9, 12, 15}, {19, 26, 33}}));
    CHECK(AC.transpose() == mnp::Matrix({{9, 19}, {12, 26}, {15, 33}}));
    CHECK(AC.transpose() == C.transpose() * A.transpose());
  }

  SUBCASE("4. A × (1, 1)") {
    const std::vector<double> un_un = {1, 1};
    CHECK(A * un_un == mnp::Matrix({{3}, {7}}));
  }
}

TEST_CASE("Matrice : (AB)^T = B^T A^T sur des matrices non carrées") {
  SUBCASE("(2x3)(3x4)") {
    const mnp::Matrix a({{1, -2, 3}, {0.5, 4, -1}});
    const mnp::Matrix b({{2, 0, 1, -3}, {1, 1, 0, 2}, {-1, 5, 2, 0.25}});
    const mnp::Matrix abt = (a * b).transpose();
    CHECK(abt.nbLignes() == 4);
    CHECK(abt.nbColonne() == 2);
    CHECK(abt == b.transpose() * a.transpose());
    CHECK_MESSAGE(a.transpose().transpose() == a, "(A^T)^T doit valoir A");
  }

  SUBCASE("(4x5)(5x3), valeurs variées") {
    mnp::Matrix a(4, 5);
    mnp::Matrix b(5, 3);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 5; j++)
        a(i, j) = (i * 7 + j * 3) % 9 - 4.5;
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 3; j++)
        b(i, j) = (i * 2 + j * 5) % 7 - 3;
    CHECK((a * b).transpose() == b.transpose() * a.transpose());
  }
}

TEST_CASE("Matrice : cas limites 1x1 et 1xn") {
  SUBCASE("1x1") {
    const mnp::Matrix trois({{3}});
    CHECK(trois * mnp::Matrix({{4}}) == mnp::Matrix({{12}}));
    CHECK(trois.transpose() == trois);
    CHECK(mnp::Matrix::identity(1) == mnp::Matrix({{1}}));
    CHECK(mnp::Matrix::identity(1).isIndentity());
    CHECK(trois * mnp::Matrix::identity(1) == trois);
    CHECK(trois * std::vector<double>{2} == mnp::Matrix({{6}}));
  }

  SUBCASE("1xn : ligne, transposée, produits scalaire et extérieur") {
    const mnp::Matrix ligne({{1, 2, 3}});
    const mnp::Matrix autre({{4, 5, 6}});
    CHECK(ligne.transpose() == mnp::Matrix({{1}, {2}, {3}}));
    // (1x3)(3x1) -> 1x1 : produit scalaire
    CHECK(ligne * autre.transpose() == mnp::Matrix({{32}}));
    // (3x1)(1x3) -> 3x3 : produit extérieur
    CHECK(ligne.transpose() * autre ==
          mnp::Matrix({{4, 5, 6}, {8, 10, 12}, {12, 15, 18}}));
    CHECK_THROWS_AS(ligne * autre, std::invalid_argument);
  }
}

TEST_CASE("Systemes lineaires") {
  const mnp::Matrix systeme({{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}});
  const std::vector<double> vector_equal = {8, -11, -3};
  std::vector<double> v = {2, 3, -1};
  const auto x = mnp::solve(systeme, vector_equal);
  REQUIRE(x.size() == 3);
  CHECK(x[0] == doctest::Approx(2));
  CHECK(x[1] == doctest::Approx(3));
  CHECK(x[2] == doctest::Approx(-1));
}
