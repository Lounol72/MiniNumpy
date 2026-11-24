// Banc d'essai de mnp::solve (Gauss) et de mnp::LuDecomposition (étape 3).
//
// Usage : ./bench_solve [-k K] [n1 n2 ...]   (par défaut : K = 100, n = 200 500 1000)
// À lancer depuis le build Release, sinon les sanitizers faussent les temps :
//   cmake -B build-release -DCMAKE_BUILD_TYPE=Release
//   cmake --build build-release && ./build-release/bench_solve
//
// Pour chaque n : une matrice aléatoire n × n et un second membre b, tirés
// avec une graine fixe (mêmes données à chaque lancement, donc mesures
// comparables d'une version de solve à l'autre). solve est lancé plusieurs
// fois et on garde le temps le plus court : c'est le moins perturbé par le
// reste du système.
//
// Premier tableau : un seul b. Gauss (solve) contre LU (décomposition +
// resoudre) : les deux font ~(2/3)·n³ opérations, les temps doivent être proches.
//
// Deuxième tableau : la même A avec K seconds membres différents (K = 100 par
// défaut). C'est le cas où LU doit gagner : solve refait toute l'élimination
// (~(2/3)·n³) pour chaque b, alors que LU décompose A une seule fois puis ne
// fait qu'une descente + remontée (~2n²) par b. Les résolutions sont mesurées
// de deux façons : K appels à resoudre(b), qui relisent L et U K fois, et un
// seul resoudreMatrice(B), qui lit chaque case de L et U une seule fois.
// Temps en ms.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "mnp/Linalg.hpp"
#include "mnp/LuDecomposition.hpp"
#include "mnp/Matrix.hpp"

namespace {

constexpr unsigned GRAINE = 42;

/// Nombre de répétitions : plus la matrice est petite, plus on répète
int repetitions(int n) {
  if (n <= 200)
    return 20;
  if (n <= 500)
    return 5;
  return 3;
}

struct Systeme {
  mnp::Matrix A;
  std::vector<double> b;
};

/// Matrice et second membre aléatoires, coefficients uniformes dans [-1, 1]
Systeme systeme_aleatoire(int n, std::mt19937 &rng) {
  std::uniform_real_distribution<double> loi(-1.0, 1.0);
  Systeme s{mnp::Matrix(n, n), std::vector<double>(n)};
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      s.A(i, j) = loi(rng);
    s.b[i] = loi(rng);
  }
  return s;
}

/// K seconds membres aléatoires, à la suite de A dans la même suite aléatoire
std::vector<std::vector<double>> seconds_membres(int n, int k,
                                                 std::mt19937 &rng) {
  std::uniform_real_distribution<double> loi(-1.0, 1.0);
  std::vector<std::vector<double>> bs(k, std::vector<double>(n));
  for (auto &b : bs)
    for (double &v : b)
      v = loi(rng);
  return bs;
}

/// Nombre de répétitions pour K résolutions : une seule au-delà de n = 200,
/// sinon n = 1000 prendrait plusieurs dizaines de secondes
int repetitions_k(int n) { return n <= 200 ? 3 : 1; }

/// Durée en millisecondes depuis debut
double ms_depuis(std::chrono::steady_clock::time_point debut) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - debut)
      .count();
}

/// max |(Ax − b)ᵢ| : 0 pour une solution exacte, ~1e-12 est très bon
double residu_max(const mnp::Matrix &A, const std::vector<double> &x,
                  const std::vector<double> &b) {
  const mnp::Matrix Ax = A * x; // matrice n × 1 (décision de l'étape 2)
  double max = 0.0;
  for (int i = 0; i < A.nbLignes(); i++)
    max = std::max(max, std::fabs(Ax(i, 0) - b[i]));
  return max;
}

} // namespace

int main(int argc, char *argv[]) {
#ifndef NDEBUG
  std::cerr
      << "Attention : build Debug (assert et sanitizers actifs), les "
         "temps ne sont pas représentatifs. Utiliser le build Release.\n\n";
#endif

  int K = 100;
  std::vector<int> tailles;
  for (int k = 1; k < argc; k++) {
    if (std::string(argv[k]) == "-k" && k + 1 < argc)
      K = std::atoi(argv[++k]);
    else
      tailles.push_back(std::atoi(argv[k]));
  }
  if (K <= 0) {
    std::cerr << "K doit être positif\n";
    return 1;
  }
  if (tailles.empty())
    tailles = {200, 500, 1000};

  std::cout << "Un seul second membre\n"
            << std::setw(6) << "n" << std::setw(14) << "Gauss (ms)"
            << std::setw(12) << "GFLOP/s" << std::setw(14) << "résidu"
            << std::setw(12) << "LU (ms)" << std::setw(12) << "GFLOP/s"
            << std::setw(14) << "résidu" << "\n";

  for (int n : tailles) {
    if (n <= 0) {
      std::cerr << "Taille ignorée : " << n << "\n";
      continue;
    }
    std::mt19937 rng(GRAINE); // même graine pour chaque n : reproductible
    const Systeme s = systeme_aleatoire(n, rng);

    double gauss_ms = 0.0, lu_ms = 0.0;
    std::vector<double> x_gauss, x_lu;
    for (int r = 0; r < repetitions(n); r++) {
      auto debut = std::chrono::steady_clock::now();
      x_gauss = mnp::solve(s.A, s.b);
      const double ms_g = ms_depuis(debut);

      debut = std::chrono::steady_clock::now();
      x_lu = mnp::LuDecomposition(s.A).resoudre(s.b);
      const double ms_l = ms_depuis(debut);

      if (r == 0 || ms_g < gauss_ms)
        gauss_ms = ms_g;
      if (r == 0 || ms_l < lu_ms)
        lu_ms = ms_l;
    }

    // Les deux méthodes coûtent environ (2/3)·n³ opérations
    const double operations = 2.0 / 3.0 * n * static_cast<double>(n) * n;
    auto gflops = [&](double ms) { return operations / (ms * 1e-3) / 1e9; };

    std::cout << std::setw(6) << n << std::fixed << std::setprecision(2)
              << std::setw(14) << gauss_ms << std::setw(12) << gflops(gauss_ms)
              << std::scientific << std::setprecision(1) << std::setw(14)
              << residu_max(s.A, x_gauss, s.b) << std::fixed
              << std::setprecision(2) << std::setw(12) << lu_ms
              << std::setw(12) << gflops(lu_ms) << std::scientific
              << std::setprecision(1) << std::setw(14)
              << residu_max(s.A, x_lu, s.b) << std::defaultfloat << "\n";
  }

  // ==========================================================================
  // Même A, K seconds membres
  // ==========================================================================
  std::cout << "\nMême A, " << K << " seconds membres différents"
            << " (temps de résolution sans la décomposition)\n"
            << std::setw(6) << "n" << std::setw(12) << "Gauss" << std::setw(12)
            << "décomp LU" << std::setw(14) << "K×resoudre" << std::setw(18)
            << "resoudreMatrice" << std::setw(10) << "gain" << std::setw(12)
            << "LU total" << std::setw(12) << "vs Gauss" << std::setw(12)
            << "résidu G" << std::setw(12) << "résidu b" << std::setw(12)
            << "résidu B" << "\n";

  for (int n : tailles) {
    if (n <= 0)
      continue;
    std::mt19937 rng(GRAINE);
    const Systeme s = systeme_aleatoire(n, rng); // même A que plus haut
    const auto bs = seconds_membres(n, K, rng);

    // Les mêmes K seconds membres rangés en colonnes : B est n × K
    mnp::Matrix B(n, K);
    for (int k = 0; k < K; k++)
      for (int i = 0; i < n; i++)
        B(i, k) = bs[k][i];

    double gauss_ms = 0.0, decomp_ms = 0.0, vect_ms = 0.0, mat_ms = 0.0;
    std::vector<std::vector<double>> xs_gauss(K), xs_lu(K);
    mnp::Matrix X(n, K);
    for (int r = 0; r < repetitions_k(n); r++) {
      // Gauss : une élimination complète par b
      auto debut = std::chrono::steady_clock::now();
      for (int k = 0; k < K; k++)
        xs_gauss[k] = mnp::solve(s.A, bs[k]);
      const double ms_g = ms_depuis(debut);

      // LU : A décomposée une seule fois
      debut = std::chrono::steady_clock::now();
      const mnp::LuDecomposition lu(s.A);
      const double ms_d = ms_depuis(debut);

      // puis K descentes + remontées, un b à la fois : L et U relues K fois
      debut = std::chrono::steady_clock::now();
      for (int k = 0; k < K; k++)
        xs_lu[k] = lu.resoudre(bs[k]);
      const double ms_v = ms_depuis(debut);

      // ou tous les b d'un coup : chaque case de L et U lue une seule fois
      debut = std::chrono::steady_clock::now();
      X = lu.resoudreMatrice(B);
      const double ms_m = ms_depuis(debut);

      if (r == 0 || ms_g < gauss_ms)
        gauss_ms = ms_g;
      if (r == 0 || ms_d < decomp_ms)
        decomp_ms = ms_d;
      if (r == 0 || ms_v < vect_ms)
        vect_ms = ms_v;
      if (r == 0 || ms_m < mat_ms)
        mat_ms = ms_m;
    }

    double residu_gauss = 0.0, residu_b = 0.0, residu_B = 0.0;
    for (int k = 0; k < K; k++) {
      residu_gauss = std::max(residu_gauss, residu_max(s.A, xs_gauss[k], bs[k]));
      residu_b = std::max(residu_b, residu_max(s.A, xs_lu[k], bs[k]));
      std::vector<double> xk(n);
      for (int i = 0; i < n; i++)
        xk[i] = X(i, k);
      residu_B = std::max(residu_B, residu_max(s.A, xk, bs[k]));
    }

    const double lu_total = decomp_ms + mat_ms;
    std::cout << std::setw(6) << n << std::fixed << std::setprecision(1)
              << std::setw(12) << gauss_ms << std::setw(12) << decomp_ms
              << std::setw(14) << vect_ms << std::setw(18) << mat_ms
              << std::setw(9) << vect_ms / mat_ms << "×" << std::setw(12)
              << lu_total << std::setw(11) << gauss_ms / lu_total << "×"
              << std::scientific << std::setprecision(1) << std::setw(12)
              << residu_gauss << std::setw(12) << residu_b << std::setw(12)
              << residu_B << std::defaultfloat << "\n";
  }

  return 0;
}
