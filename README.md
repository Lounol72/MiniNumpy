# mini-numpy

Bibliothèque de calcul numérique écrite de zéro en C++20, pour réapprendre le C++ et les maths du Master.
Le but final : réimplémenter en C++, avec cette bibliothèque à la place de NumPy, le réseau de neurones du livre *Neural Networks from Scratch in Python*, jusqu'à l'entraîner sur Fashion-MNIST.

- `include/mnp/` : en-têtes publics (`#include "mnp/..."`)
- `src/` : implémentation de la bibliothèque
- `apps/` : programmes de démonstration (un par étape)
- `tests/` : tests doctest
- `data/` : jeux de données (CSV, Fashion-MNIST)
- `third_party/` : dépendances embarquées (doctest)

Feuille de route et journal : note Obsidian `Dev/00 - Mini-NumPy.md`.

## Référence

Le projet final s'appuie sur le livre suivant :

> Harrison Kinsley et Daniel Kukieła, *Neural Networks from Scratch in Python*, 2020. <https://nnfs.io>

Le code de ce dépôt est une réécriture personnelle en C++. Il ne reprend ni le texte ni le code Python du livre.
