// Banc d'essai de mnp::solve (étape 3, section 7 du guide).
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
// Deuxième tableau : la même A avec K seconds membres différents (K = 100 par
// défaut). C'est le cas où LU doit gagner : aujourd'hui solve refait toute
// l'élimination (~(2/3)·n³) pour chaque b, alors que seul b change.
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

  std::cout << std::setw(6) << "n" << std::setw(14) << "temps (ms)"
            << std::setw(12) << "GFLOP/s" << std::setw(16) << "résidu max"
            << "\n";

  for (int n : tailles) {
    if (n <= 0) {
      std::cerr << "Taille ignorée : " << n << "\n";
      continue;
    }
    std::mt19937 rng(GRAINE); // même graine pour chaque n : reproductible
    const Systeme s = systeme_aleatoire(n, rng);

    double meilleur_ms = 0.0;
    std::vector<double> x;
    for (int r = 0; r < repetitions(n); r++) {
      const auto debut = std::chrono::steady_clock::now();
      x = mnp::solve(s.A, s.b);
      const auto fin = std::chrono::steady_clock::now();
      const double ms =
          std::chrono::duration<double, std::milli>(fin - debut).count();
      if (r == 0 || ms < meilleur_ms)
        meilleur_ms = ms;
    }

    // L'élimination de Gauss coûte environ (2/3)·n³ opérations
    const double operations = 2.0 / 3.0 * n * static_cast<double>(n) * n;
    const double gflops = operations / (meilleur_ms * 1e-3) / 1e9;

    std::cout << std::setw(6) << n << std::setw(14) << std::fixed
              << std::setprecision(2) << meilleur_ms << std::setw(12)
              << std::setprecision(2) << gflops << std::setw(16)
              << std::scientific << std::setprecision(1)
              << residu_max(s.A, x, s.b) << std::defaultfloat << "\n";
  }

  // ==========================================================================
  // Même A, K seconds membres
  // ==========================================================================
  std::cout << "\nMême A, " << K << " seconds membres différents\n"
            << std::setw(6) << "n" << std::setw(16) << "total (ms)"
            << std::setw(14) << "par b (ms)" << std::setw(16) << "résidu max"
            << "\n";

  for (int n : tailles) {
    if (n <= 0)
      continue;
    std::mt19937 rng(GRAINE);
    const Systeme s = systeme_aleatoire(n, rng); // même A que plus haut
    const auto bs = seconds_membres(n, K, rng);

    double meilleur_ms = 0.0;
    std::vector<std::vector<double>> xs(K);
    for (int r = 0; r < repetitions_k(n); r++) {
      const auto debut = std::chrono::steady_clock::now();
      // Aujourd'hui : une élimination complète par b.
      // Avec LU : décomposer A une seule fois ici, avant la boucle, puis
      // descente + remontée pour chaque b.
      for (int j = 0; j < K; j++)
        xs[j] = mnp::solve(s.A, bs[j]);
      const auto fin = std::chrono::steady_clock::now();
      const double ms =
          std::chrono::duration<double, std::milli>(fin - debut).count();
      if (r == 0 || ms < meilleur_ms)
        meilleur_ms = ms;
    }

    double residu = 0.0;
    for (int j = 0; j < K; j++)
      residu = std::max(residu, residu_max(s.A, xs[j], bs[j]));

    std::cout << std::setw(6) << n << std::setw(16) << std::fixed
              << std::setprecision(1) << meilleur_ms << std::setw(14)
              << std::setprecision(2) << meilleur_ms / K << std::setw(16)
              << std::scientific << std::setprecision(1) << residu
              << std::defaultfloat << "\n";
  }

  return 0;
}
