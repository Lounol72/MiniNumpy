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

  // Seul le produit est chronométré, pas le remplissage.
  // (Comparaison i,j,k / i,k,j du 16/11 : voir le journal de bord.)
  const auto debut = std::chrono::steady_clock::now();
  const mnp::Matrix c = a * b;
  const auto fin = std::chrono::steady_clock::now();

  const double ms =
      std::chrono::duration<double, std::milli>(fin - debut).count();
  const double operations = 2.0 * n * n * n;
  std::cout << "Produit " << n << "x" << n << " : " << ms << " ms, soit "
            << operations / (ms / 1000.0) / 1e6
            << " millions d'opérations par seconde\n";

  // Utiliser le résultat (avec -O2, le compilateur pourrait sinon supprimer
  // un calcul dont personne ne lit le résultat)
  std::cout << "c(1, 1) = " << c(1, 1) << '\n';
  return 0;
}
