#include <chrono>
#include <iostream>

#include "mnp/Matrix.hpp"

int main(void) {
  // Affichage
  std::vector<std::vector<double>> v = {{1, 2, 3}, {4, 5, 6}};
  std::vector<std::vector<double>> y = {{1, -22.5}, {333, 0}};
  mnp::Matrix mata(v);
  mnp::Matrix matb(y);
  std::cout << mata << '\n';
  std::cout << matb << "\n\n";

  // Mesure du produit n×n
  const int n = 500;
  mnp::Matrix a(n, n);
  mnp::Matrix b(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      a(i, j) = (i + j) % 10;
      b(i, j) = (i * j) % 7;
    }

  // Seul le produit est chronométré, pas le remplissage
  const double operations = 2.0 * n * n * n;
  auto mesurer = [&](const char *nom, auto produit) {
    const auto debut = std::chrono::steady_clock::now();
    mnp::Matrix c = produit();
    const auto fin = std::chrono::steady_clock::now();
    const double ms =
        std::chrono::duration<double, std::milli>(fin - debut).count();
    std::cout << "Produit " << n << "x" << n << " (ordre " << nom
              << ") : " << ms << " ms, soit "
              << operations / (ms / 1000.0) / 1e6
              << " millions d'opérations par seconde\n";
    return c;
  };

  const mnp::Matrix c_ijk = mesurer("i,j,k", [&] { return a.produit_ijk(b); });
  const mnp::Matrix c_ikj = mesurer("i,k,j", [&] { return a * b; });

  // Utiliser les résultats (à l'étape 6, avec -O2, le compilateur pourrait
  // sinon supprimer un calcul dont personne ne lit le résultat)
  std::cout << "Mêmes résultats : " << (c_ijk == c_ikj ? "oui" : "NON") << '\n';
  return 0;
}
