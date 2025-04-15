# BE Elasticité

Ce dépôt permet de réaliser des études d’élasticité 2D par éléments finis en Python, avec une structure modulaire pour comparer facilement différents paramètres physiques (comme le coefficient de Poisson) ou numériques (raffinement du maillage).

Le cœur du projet est la classe `ElasticitySolver` dans elasticite.py, refactorisée à partir de l’ancien script main_elasticite_stat_2D_propre.py. Cette classe peut être utilisée directement dans le notebook elasticite.ipynb pour lancer des calculs, visualiser les résultats et comparer différentes configurations.

---

## Installation et environnement

Ce projet utilise [Poetry](https://python-poetry.org/) pour la gestion des dépendances et de l’environnement virtuel.

1. **Cloner le dépôt :**
   ```bash
   git clone https://github.com/votre-utilisateur/BE_Elasticite.git
   cd BE_Elasticite
   ```

2. **Installer Poetry (si besoin) :**

> https://python-poetry.org/docs/#installing-with-the-official-installer

> ***NB** : Ne pas oublier d'ajouter poetry au PATH.*

3. **Créer l’environnement et installer les dépendances :**
   ```bash
   poetry install
   ```

   > **Remarque :**  
   La version de `matplotlib` doit est **3.6** pour garantir la compatibilité avec le code.

4. **Activer l’environnement Poetry :**
   ```bash
   poetry env activate
   ```

5. **Optionnel :** Creer un Jupyter Kernel
    ```bash
    poetry run python -m ipykernel install --user --name=Be-Elasticite --display-name "Python (Be-Elasticite)"
    ````
---

## Utilisation

- **Classe principale :**  
  La classe `ElasticitySolver` (dans elasticite.py) permet de :
  - Charger un maillage (choix parmi plusieurs raffinements ou géométries)
  - Définir les propriétés mécaniques (module d’Young, coefficient de Poisson)
  - Résoudre le problème d’élasticité 2D
  - Visualiser le maillage initial, déformé et les contraintes de Von Mises

- **Exemple dans un notebook (`elasticite.ipynb`) :**
  ```python
  from elasticite import ElasticitySolver

  # Création et résolution pour un maillage et un coefficient de Poisson donnés
  solver = ElasticitySolver(maille=2, v=0.3)
  solver.solve()
  solver.plot_results()

  # Comparaison pour différents coefficients de Poisson
  for v in [0.0, 0.3, 0.49]:
      solver = ElasticitySolver(maille=2, v=v)
      solver.solve()
      # Accès aux résultats : solver.results['U'], solver.results['SVM'], etc.
  ```

- **Comparaison des paramètres :**
  - Modifiez `maille` pour changer le raffinement ou la géométrie.
  - Modifiez `v` pour comparer l’effet du coefficient de Poisson.

---

## Structure du projet

```
BE_Elasticite/
├── elasticite.py                # Classe principale (modulaire)
├── elasticite.ipynb             # Notebook d’utilisation et de comparaison
├── main_elasticite_stat_2D_propre.py # Ancien script (pour référence)
├── pyproject.toml               # Dépendances Poetry
├── maillage_*.mat               # Fichiers de maillage (nécessaires)
└── README.md
```



---

## Licence

Projet académique – usage pédagogique uniquement.

---

N’hésitez pas à ouvrir une issue pour toute question ou suggestion !