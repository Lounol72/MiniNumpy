#pragma once

#include "mnp/Matrix.hpp"
#include "mnp/Precision.hpp"
#include <vector>

namespace mnp {

/**
 * @brief Trouve la solution pour la matrice m selon le vecteur v
 *
 * @param m mnp::Matrix
 * @param v std::vector
 */
std::vector<double> solve(const mnp::Matrix &m, const std::vector<double> &v);

/**
 * @brief Pivot de Gauss sur matrice m
 *
 * @param m mnp::Matrix
 */
mnp::Matrix pivotDeGauss(const mnp::Matrix &m);

/**
 * @brief Algorithme de remontée permettant de résoudre un système d'équations linéaire
 *
 * @param m Matrice triangulaire
 */
std::vector<double> remontee(const mnp::Matrix &m);

} // namespace mnp
