# Simulation de Ponts Thermiques

Ce projet permet de simuler et d'analyser les ponts thermiques dans différentes configurations de murs et planchers. Il utilise la méthode des différences finies et le schéma de Jacobi pour résoudre l'équation de la chaleur en 2D.

## Description

Le projet permet d'étudier différentes configurations constructives pour réduire les pertes thermiques au niveau des jonctions mur/plancher:

- Mur simple avec isolation 
- Mur avec rupteur thermique
- Mur avec planelle isolante
- Différentes positions d'isolation

## Fonctionnalités

- Simulation thermique 2D en régime stationnaire
- Calcul des champs de température
- Visualisation des flux thermiques
- Calcul du coefficient de pont thermique Ψ
- Estimation des pertes énergétiques
- Comparaison de différentes solutions constructives

## Structure du code

- `main.py` : Programme principal et paramètres de simulation
- `geometrie.py` : Définition des différentes configurations géométriques
- `fonction_calcul.py` : Fonctions de calcul numérique (Jacobi, gradients...)
- `fonction_therm.py` : Fonctions thermiques (conditions limites, flux...)

## Paramètres ajustables

- Dimensions géométriques (épaisseurs, hauteurs...)
- Propriétés des matériaux (conductivités thermiques...)
- Conditions aux limites (températures intérieure/extérieure)
- Paramètres numériques (discrétisation, relaxation...)

## Utilisation

1. Définir les paramètres dans `main.py`
2. Choisir une configuration géométrique 
3. Lancer la simulation
4. Visualiser les résultats (champs de température, flux)
5. Analyser les performances (coefficient Ψ, coût énergétique)

## Résultats

Le programme génère:
- Une visualisation de la géométrie
![Géométrie du pont thermique]
<img src="planelle.png" width="300"/>
- Les champs de température
- Les vecteurs de flux thermique
![Température et flux thermique]
<img src="temp_planelle.png" width="300"/>
- Le coefficient de pont thermique Ψ
- Une estimation des pertes énergétiques

## Dépendances

- NumPy: Calculs numériques
- Matplotlib: Visualisation des résultats

