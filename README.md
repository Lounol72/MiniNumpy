# MiniNumpy

Bibliothèque de calcul numérique écrite de zéro en C++20, pour réapprendre le C++ et les maths du Master.

Le but final : réimplémenter en C++ le réseau de neurones du livre *Neural Networks from Scratch in Python*, avec cette bibliothèque à la place de NumPy, jusqu'à l'entraîner sur Fashion-MNIST.

Le projet avance par étapes. Chaque étape ajoute une notion de maths, une notion de C++ et un outil de build (Make, puis CMake).

## Avancement

| Étape | Sujet | État |
|---|---|---|
| 0 | Compiler à la main, Makefile, doctest | ✅ |
| 1 | `Vector` : norme, produit scalaire, angle | ✅ |
| 2 | `Matrix` : opérations, produit, transposée, affichage | ✅ (relecture en cours) |
| 3 | Systèmes linéaires : Gauss avec pivot partiel, décomposition LU | 🚧 en cours |
| 4 | Moindres carrés et régression | à venir |
| 5 | Aléatoire et Monte Carlo | à venir |
| 6 | Optimisation : descente de gradient, Newton | à venir |
| 7 | Valeurs propres et ACP | à venir |
| 8 | Différentiation automatique | à venir |
| 9 | Projet final : le réseau de neurones de NNFS sur Fashion-MNIST | à venir |

À l'étape 3, ce qui marche déjà :
- `mnp::solve` : élimination de Gauss avec pivot partiel
- `mnp::LuDecomposition` : décomposition PA = LU, `resoudre(b)` et `resoudreMatrice(B)` pour plusieurs seconds membres

Il reste à écrire `determinant()` et `inverse()` à partir de LU.

## Compiler

Prérequis : un compilateur C++20 (g++ ou clang++) et **CMake ≥ 4.0**.

```sh
# Debug (par défaut) : -Werror, sanitizers address et undefined
cmake -B build
cmake --build build
ctest --test-dir build

# Release : pour les mesures de temps, sans sanitizers
cmake -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

CMake compile chaque fichier de `apps/` en un exécutable, et chaque fichier de `tests/` en un programme de test lancé par `ctest`. Le fichier `compile_commands.json` (pour clangd) est généré dans `build/`.

Le `Makefile` vient des étapes 0 à 2. Il est gardé le temps de finir la migration vers CMake. Il sert encore à générer le carnet de maths (`make doc`).

## Exemple

Résoudre le système suivant :

```
 2x +  y −  z =   8
−3x −  y + 2z = −11
−2x +  y + 2z =  −3
```

```cpp
#include "mnp/Linalg.hpp"
#include "mnp/LuDecomposition.hpp"

const mnp::Matrix A({{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}});
const std::vector<double> b = {8, -11, -3};

// Gauss : une élimination à chaque appel
auto x = mnp::solve(A, b);             // {2, 3, -1}

// LU : on décompose A une seule fois, puis chaque résolution coûte ~2n² opérations
const mnp::LuDecomposition lu(A);
auto y = lu.resoudre(b);               // {2, 3, -1}
```

## Démos

Après `cmake --build build`, les exécutables sont dans `build/` :

| Programme | Ce qu'il montre |
|---|---|
| `main` | la version de la bibliothèque (étape 0) |
| `main_vector` | l'angle entre deux vecteurs (étape 1) |
| `main_matrix` | l'affichage des matrices et le temps d'un produit n × n (étape 2) |
| `main_linalg` | la résolution d'un système 3 × 3 par Gauss (étape 3) |
| `bench_solve` | Gauss contre LU, avec un ou plusieurs seconds membres (étape 3) |

`bench_solve` se lance depuis le build Release : `./build-release/bench_solve [-k K] [n1 n2 ...]`.

## Organisation

```
include/mnp/   en-têtes publics (#include "mnp/...")
src/           implémentation de la bibliothèque (libmnp)
apps/          programmes de démonstration
tests/         tests doctest
docs/carnet/   carnet de maths en LaTeX
data/          jeux de données (CSV, Fashion-MNIST)
third_party/   dépendances embarquées (doctest)
```

## Choix de conception

- **Stockage** : une `Matrix` range ses valeurs ligne par ligne dans un seul `std::vector<double>`. La case (i, j) est à l'indice `i * colonnes + j`.
- **Comparaisons de `double`** : jamais avec `==`. `mnp::approx_equal` (dans `Precision.hpp`) compare avec une tolérance relative, `PRECISION = 1e-9` par défaut.
- **Accès** : `m(i, j)` ne vérifie rien, `m.at(i, j)` vérifie les indices et lève une exception.
- **Fabriques statiques** : `Matrix::zeros`, `Matrix::ones`, `Matrix::identity`, `Vector::zeros`, `Vector::ones`.
- **Itérateurs** : `begin()` et `end()` permettent les boucles `for` sur intervalle et les algorithmes de `<algorithm>`.
- **Rigueur** : `-Wall -Wextra -Wpedantic -Wshadow -Werror` partout, sanitizers en Debug. Chaque fonction a un test écrit avant le code.
- **Langue** : les noms sont en français (`nbLignes`, `resoudre`, `ajouterColonne`), sauf quand ils reprennent un nom connu de NumPy (`zeros`, `ones`, `transpose`).

## Limites

- L'API ne cherche pas à imiter celle de NumPy : pas de broadcasting, pas de tableaux à plus de deux dimensions.

## Carnet de maths

`docs/carnet/` contient un carnet en LaTeX, avec un chapitre par étape : les maths derrière chaque algorithme, faites à la main sur des exemples. Pour générer le PDF (latexmk nécessaire) :

```sh
make doc    # produit docs/carnet/carnet_maths.pdf
```

## Référence

Le projet final s'appuie sur le livre suivant :

> Harrison Kinsley et Daniel Kukieła, *Neural Networks from Scratch in Python*, 2020. <https://nnfs.io>

Le code de ce dépôt est une réécriture personnelle en C++. Il ne reprend ni le texte ni le code Python du livre.
