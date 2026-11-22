#pragma once

#include <ostream>
#include <vector>

#include "mnp/Precision.hpp"
#include "mnp/Vector.hpp"

namespace mnp {

class Matrix {
public:
  // ========================================================================
  // Constructeurs
  // ========================================================================

  /**
   * @brief Crée une matrice i × j remplie de 0
   *
   * @param i nombre de lignes
   * @param j nombre de colonnes
   */
  explicit Matrix(int i, int j);

  /**
   * @brief Crée une matrice à partir d'une liste de lignes
   *
   * @param vect les lignes, toutes de la même longueur
   * @throw std::invalid_argument si vect est vide, si une ligne est vide ou
   * si les lignes n'ont pas toutes la même longueur
   */
  explicit Matrix(const std::vector<std::vector<double>> &vect);

  /**
   * @brief Matrice identité n × n
   */
  static Matrix identity(int n);

  // ========================================================================
  // Propriétés
  // ========================================================================

  /// Vrai si la matrice est carrée (n × n)
  bool isSquare() const;
  /// Vrai si tous les éléments sont nuls (à PRECISION près)
  bool isZero() const;
  /// Vrai si la matrice est l'identité (à PRECISION près)
  bool isIdentity() const;

  /// Transposée : (i, j) → (j, i)
  Matrix transpose() const;

  Matrix &ajouterColonne(const mnp::Vector v);
  Matrix &ajouterLigne(const mnp::Vector v);
  int nbLignes() const { return ligne; }
  int nbColonne() const { return colonne; }
  std::vector<double> getLigne(int index) const;
  std::vector<double> getColonne(int index) const;

  // ========================================================================
  // Accès aux éléments
  // ========================================================================

  /**
   * @brief Élément (i, j)
   * @return double
   */
  double &operator()(int i, int j) { return m_values[idx(i, j)]; }
  const double &operator()(int i, int j) const { return m_values[idx(i, j)]; }

  /**
   * @brief Élément (i, j)
   * @throw std::invalid_argument si (i, j) est hors de la matrice
   */
  double &at(int i, int j);
  const double &at(int i, int j) const;

  // ========================================================================
  // Opérateurs
  // ========================================================================

  Matrix &operator+=(const Matrix &other);
  Matrix &operator-=(const Matrix &other);
  friend Matrix operator+(Matrix lhs, const Matrix &rhs);
  friend Matrix operator-(Matrix lhs, const Matrix &rhs);

  /**
   * @brief Produit matriciel, ordre de boucles i, k, j
   *
   * Pour chaque ligne i de A et chaque k, A(i,k) est lu une seule fois puis
   * multiplié par toute la ligne k de B, ajoutée à la ligne i du résultat.
   * Dans la boucle intérieure, B et le résultat sont parcourus case par case,
   * dans l'ordre où ils sont rangés en mémoire.
   *
   * @throw std::invalid_argument si colonnes de A != lignes de B
   */
  Matrix operator*(const Matrix &other) const;

  /**
   * @brief Matrice × vecteur : le vecteur est vu comme une colonne
   *
   * @return une Matrix n × 1 (choix de conception : tout est matrice)
   * @throw std::invalid_argument si la taille du vecteur != colonnes de m
   */
  friend Matrix operator*(const Matrix &m, const std::vector<double> &vect);

  Matrix operator+() const;
  Matrix operator-() const;

  /**
   * @brief Égalité élément par élément, à PRECISION près
   * @throw std::invalid_argument si les tailles diffèrent
   */
  bool operator==(const Matrix &m) const;

  /**
   * @brief Affiche la matrice à la manière de NumPy, colonnes alignées
   *
   * L'état du flux (précision, flags) est restauré après l'affichage.
   */
  friend std::ostream &operator<<(std::ostream &out, const Matrix &m);

  // ========================================================================
  // Itérateurs
  // ========================================================================

  using iterator = std::vector<double>::iterator;
  using const_iterator = std::vector<double>::const_iterator;

  /**
   * @brief Parcours de tous les éléments, dans l'ordre du stockage par lignes
   *
   * (0,0), (0,1), …, (0,n-1), puis (1,0), … : pour [[1, 2, 3], [4, 5, 6]], la
   * boucle `for (double x : m)` lit 1, 2, 3, 4, 5, 6. Ce n'est pas le choix de
   * NumPy, où `for ligne in A` parcourt les lignes. Cet ordre fait partie du
   * contrat (il est testé) : changer de stockage le changerait.
   *
   * Les itérateurs restent valides tant que la taille de la matrice ne change
   * pas. Une affectation `m = autre` les invalide (le tableau est remplacé),
   * de même que tout ajout de ligne ou de colonne.
   */
  iterator begin() { return m_values.begin(); }
  iterator end() { return m_values.end(); }
  const_iterator begin() const { return m_values.begin(); }
  const_iterator end() const { return m_values.end(); }

private:
  /// Position de l'élément (i, j) dans m_values (rangement par lignes)
  int idx(int i, int j) const { return i * colonne + j; }

  int ligne;
  int colonne;
  std::vector<double> m_values;
};

} // namespace mnp
