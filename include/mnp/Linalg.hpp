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
 * @brief Élimination de Gauss avec pivot partiel, sur place
 *
 * À chaque colonne i, choisit comme pivot l'élément de plus grande valeur
 * absolue des lignes i à n-1, échange les lignes, puis annule les éléments
 * sous le pivot. À la fin, m est triangulaire supérieure (sur ses
 * nbLignes() premières colonnes). Les colonnes en plus (par exemple le
 * second membre b d'une matrice augmentée) subissent les mêmes opérations.
 *
 * @param[in,out] m Matrice à au moins nbLignes() colonnes, modifiée sur place
 * @throw std::domain_error si un pivot est nul (au sens de approx_equal) :
 *        matrice singulière. m est alors laissée à moitié transformée.
 */
void pivotDeGauss(mnp::Matrix &m);

/**
 * @brief Algorithme de remontée permettant de résoudre un système d'équations linéaire
 *
 * @param m Matrice triangulaire
 */
std::vector<double> remontee(const mnp::Matrix &m);

} // namespace mnp
