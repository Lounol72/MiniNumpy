#pragma once

#include <algorithm>
#include <cmath>

namespace mnp {

/// Tolérance par défaut des comparaisons de double, partagée par Vector et
/// Matrix. Une vraie constante (typée, dans le namespace) plutôt qu'un
/// #define : pas de collision avec une PRECISION de l'utilisateur.
inline constexpr double PRECISION = 1e-9;

/**
 * @brief Vrai si a et b sont égaux à eps près, en tolérance relative
 *
 * L'écart est comparé à eps × max(1, |a|, |b|) : en absolu près de 0, en
 * relatif pour les grandes valeurs (1e12 et 1e12 + 1 sont alors égaux).
 */
inline bool approx_equal(double a, double b, double eps = PRECISION) {
  const double echelle = std::max({1.0, std::fabs(a), std::fabs(b)});
  return std::fabs(a - b) <= eps * echelle;
}

} // namespace mnp
