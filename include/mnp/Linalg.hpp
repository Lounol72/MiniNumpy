#pragma once

#include "mnp/Matrix.hpp"
#include <vector>

namespace mnp {

/**
 * @brief Trouve la solution pour la matrice m selon le vecteur v
 *
 * @param m mnp::Matrix
 * @param v std::vector
 */
std::vector<double> solve(const mnp::Matrix &m, const std::vector<double> &v);

} // namespace mnp
