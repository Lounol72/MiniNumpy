#pragma once

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>

namespace mnp {

class Vector {
public:
  // ========================================================================
  // Constructeur
  // ========================================================================

  /// Prend le std::vector par valeur puis le déplace : une seule copie au plus
  explicit Vector(std::vector<double> vec) : m_eles(std::move(vec)) {}

  // ========================================================================
  // Propriétés
  // ========================================================================

  std::size_t size() const { return m_eles.size(); }

  /// Norme euclidienne ‖u‖ = √(u·u)
  double norme() const;

  /**
   * @brief u / ‖u‖
   * @throw std::invalid_argument si u est le vecteur nul
   */
  Vector normalisation() const;

  /**
   * @brief Produit scalaire Σ uᵢ·xᵢ
   * @throw std::invalid_argument si les tailles diffèrent
   */
  double scalaire_analytique(const Vector &x) const;

  /**
   * @brief Angle θ ∈ [0, π] entre u et x, en radians
   *
   * cos θ est écrêté dans [-1, 1] avant acos : pour deux vecteurs colinéaires,
   * l'arrondi peut donner 1.0000000000000002, et acos renverrait NaN.
   * @throw std::invalid_argument si les tailles diffèrent ou si un vecteur est
   * nul (l'angle n'est pas défini)
   */
  double angle(const Vector &x) const;

  /**
   * @brief Produit scalaire par la formule ‖u‖·‖x‖·cos θ
   *
   * θ est lui-même calculé à partir du produit analytique : cette fonction
   * illustre la formule mais ne calcule rien de plus que scalaire_analytique.
   * Vaut 0 si l'un des vecteurs est nul.
   * @throw std::invalid_argument si les tailles diffèrent
   */
  double scalaire_geometrique(const Vector &x) const;

  // ========================================================================
  // Accès aux éléments
  // ========================================================================

  /// @throw std::out_of_range si index n'est pas dans [0, size())
  double &operator[](int index);
  double operator[](int index) const;

  // ========================================================================
  // Opérateurs composés (modifient le vecteur)
  // ========================================================================

  /// @throw std::invalid_argument si les tailles diffèrent
  Vector &operator+=(const Vector &x);
  /// @throw std::invalid_argument si les tailles diffèrent
  Vector &operator-=(const Vector &x);
  /// Produit élément par élément (comme `*` de NumPy)
  /// @throw std::invalid_argument si les tailles diffèrent
  Vector &operator*=(const Vector &x);
  Vector &operator*=(double a);
  Vector &operator/=(double a);

  // ========================================================================
  // Opérateurs binaires : forme canonique, construits à partir des composés
  // ========================================================================

  friend Vector operator+(Vector lhs, const Vector &rhs) { return lhs += rhs; }
  friend Vector operator-(Vector lhs, const Vector &rhs) { return lhs -= rhs; }
  friend Vector operator*(Vector lhs, const Vector &rhs) { return lhs *= rhs; }
  friend Vector operator*(Vector x, double a) { return x *= a; }
  friend Vector operator*(double a, Vector x) { return x *= a; }
  friend Vector operator/(Vector x, double a) { return x /= a; }
  /// a / v calcule a / vᵢ pour chaque élément (comme NumPy)
  friend Vector operator/(double a, Vector x);

  Vector operator+() const { return *this; }
  Vector operator-() const;

  /**
   * @brief Égalité élément par élément, à PRECISION près (voir approx_equal)
   *
   * Deux vecteurs de tailles différentes ne sont pas égaux : renvoie false,
   * comme std::vector et np.array_equal.
   */
  bool operator==(const Vector &x) const;

  /// Affiche {1, 2.5, -3}
  friend std::ostream &operator<<(std::ostream &os, const Vector &v);

  // ========================================================================
  // Modification de la taille
  // ========================================================================

  /**
   * @brief Ajoute val à la fin du vecteur (la taille augmente de 1)
   *
   * Peut réallouer le tableau interne : tout itérateur, pointeur ou référence
   * obtenu avant l'appel (par begin(), end() ou operator[]) devient invalide.
   * Ne pas l'appeler pendant une boucle `for (double x : v)` sur ce vecteur.
   */
  void ajouter(double val) { m_eles.push_back(val); }

  // ========================================================================
  // Itérateurs
  // ========================================================================

  using iterator = std::vector<double>::iterator;
  using const_iterator = std::vector<double>::const_iterator;

  /**
   * @brief Parcours des éléments du premier au dernier
   *
   * Permet `for (double x : v)` et les algorithmes standard
   * (std::accumulate, std::max_element…). Invalidés par ajouter().
   */
  iterator begin() { return m_eles.begin(); }
  iterator end() { return m_eles.end(); }
  const_iterator begin() const { return m_eles.begin(); }
  const_iterator end() const { return m_eles.end(); }

private:
  /// La seule vérification de taille de la classe
  /// @throw std::invalid_argument si x n'a pas la même taille que *this
  void verifierMemeTaille(const Vector &x) const;
  /// @throw std::out_of_range si index n'est pas dans [0, size())
  void verifierIndex(int index) const;

  std::vector<double> m_eles;
};

} // namespace mnp
