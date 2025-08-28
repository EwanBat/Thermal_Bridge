import numpy as np
from fonction_calcul import gradxy

############################################## Fonctions thermodynamiques #####################################################
def CL(geo, T, i, j):
    """
    Détermine les conditions limites pour un point (j, i) dans la grille
    en fonction de la géométrie 'geo' et des conditions physiques.

    Paramètres
    ----------
    geo : objet
        Contient la grille géométrique et les propriétés des milieux.
        Accès aux matériaux via geo.G[y, x] -> [type, conductivité]
    T : 2D array
        Tableau des températures.
    i : int
        Coordonnée x (indice de colonne).
    j : int
        Coordonnée y (indice de ligne).

    Retour
    ------
    tuple (bool, float)
        - bool : True si une condition limite est appliquée, sinon False.
        - float : Température calculée ou 0 si aucune condition.
    """

    # --- BORD GAUCHE : Température extérieure fixe ---
    if i == 0:
        Tij = geo.Text2
        return True, Tij

    # --- BORD DROIT : Conditions symétriques ou température intérieure ---
    elif i == geo.n - 1:
        if j == 0:  # Coin supérieur droit
            if geo.G[j, i][0] != 0:  # Non-air → symétrie
                Tij = geo.w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j + 1, i] + T[geo.m - 2, i]) \
                      + (1 - geo.w) * T[j, i]
            else:  # Air → température intérieure
                Tij = geo.Tint

        elif j == geo.m - 1:  # Coin inférieur droit
            if geo.G[j, i][0] != 0:  # Non-air → symétrie
                Tij = geo.w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j - 1, i] + T[1, i]) \
                      + (1 - geo.w) * T[j, i]
            else:  # Air → température intérieure
                Tij = geo.Tint

        else:  # Bord droit sans être dans un coin
            if geo.G[j, i][0] != 0:
                Tij = geo.w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j + 1, i] + T[j - 1, i]) \
                      + (1 - geo.w) * T[j, i]
            else:
                Tij = geo.Tint

        return True, Tij

    # --- BORD SUPÉRIEUR ---
    elif j == 0:
        Tij = geo.w * (1 / 4) * (T[j, i + 1] + T[j, i - 1] + T[j + 1, i] + T[geo.m - 2, i]) \
              + (1 - geo.w) * T[j, i]
        return True, Tij

    # --- BORD INFÉRIEUR ---
    elif j == geo.m - 1:
        Tij = geo.w * (1 / 4) * (T[j, i + 1] + T[j, i - 1] + T[j - 1, i] + T[1, i]) \
              + (1 - geo.w) * T[j, i]
        return True, Tij

    # --- INTERFACE ENTRE DEUX MATÉRIAUX SUR L'AXE X ---
    elif geo.G[j, i][0] != geo.G[j, i + 1][0] or geo.G[j, i][0] != geo.G[j, i - 1][0]:
        lamb1 = geo.G[j, i - 1][1]
        lamb2 = geo.G[j, i + 1][1]

        # Cas avec convection (type -1)
        if geo.G[j, i - 1][0] == -1 and geo.G[j, i][0] != geo.G[j, i - 1][0]:
            Tij = (lamb2 * T[j, i + 1] / geo.dx + T[j, i - 1] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        elif geo.G[j, i + 1][0] == -1 and geo.G[j, i][0] != geo.G[j, i + 1][0]:
            Tij = (lamb2 * T[j, i - 1] / geo.dx + T[j, i + 1] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        else:  # Interface purement conductrice
            Tij = (lamb1 * T[j, i - 1] + lamb2 * T[j, i + 1]) / (lamb2 + lamb1)

        return True, Tij

    # --- INTERFACE ENTRE DEUX MATÉRIAUX SUR L'AXE Y ---
    elif geo.G[j, i][0] != geo.G[j + 1, i][0] or geo.G[j, i][0] != geo.G[j - 1, i][0]:
        lamb1 = geo.G[j - 1, i][1]
        lamb2 = geo.G[j + 1, i][1]

        # Cas avec convection (type -1)
        if geo.G[j - 1, i][0] == -1 and geo.G[j, i][0] != geo.G[j - 1, i][0]:
            Tij = (lamb2 * T[j + 1, i] / geo.dx + T[j - 1, i] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        elif geo.G[j + 1, i][0] == -1 and geo.G[j, i][0] != geo.G[j + 1, i][0]:
            Tij = (lamb2 * T[j - 1, i] / geo.dx + T[j + 1, i] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        else:  # Interface purement conductrice
            Tij = (lamb1 * T[j - 1, i] + lamb2 * T[j + 1, i]) / (lamb2 + lamb1)

        return True, Tij

    # --- CAS PAR DÉFAUT : aucun traitement spécial, laisser l'itération ---
    else:
        return False, 0

def Jth(geo, T):
    """
    Calcule le champ de flux thermique J dans la grille,
    selon la loi de Fourier : 
        J = -λ * ∇T

    Paramètres
    ----------
    geo : objet
        Géométrie contenant les propriétés des matériaux.
        Accès via geo.G[y, x, :] -> [type, conductivité]
    T : 2D array
        Tableau des températures.

    Retour
    ------
    J : 3D numpy array
        Champ des flux thermiques de taille (m, n, 2)
        - J[:, :, 0] = composante selon y
        - J[:, :, 1] = composante selon x
    """
    # Initialisation du tableau des flux :
    # dimensions = [hauteur (m), largeur (n), direction (2)]
    J = np.zeros((geo.m, geo.n, 2))

    # Parcours de chaque point de la grille
    for i in range(geo.n):       # axe x
        for j in range(geo.m):   # axe y
            # Gradient local de température (∇T)
            Grad = gradxy(geo, T, i, j)

            # Application de la loi de Fourier
            J[j, i] = -geo.G[j, i, 1] * Grad

    return J

def Psi_lin(geo0, T0, geo, T, J):
    """
    Calcule le coefficient linéique Psi pour comparer les flux thermiques
    entre un mur "référence" et le pont thermique étudié.

    Paramètres
    ----------
    geo0 : objet
        Géométrie de référence (mur sans pont thermique).
        Accès aux matériaux via geo0.G[y, x, :] -> [type, conductivité]
    T0 : 2D array
        Tableau des températures dans le cas du mur seul.
    geo : objet
        Géométrie du modèle avec pont thermique.
    T : 2D array
        Tableau des températures dans le cas avec pont thermique.
    J : 3D array
        Tableau des flux thermiques calculés (J[y, x, direction]).
        direction = 0 pour y, 1 pour x.

    Retour
    ------
    float
        Valeur du coefficient linéique Psi (W/m·K).
    """

    # --- Initialisation ---
    i = 1  # On commence après l'air extérieur
    L_U0j = []        # Flux conductifs par segment du mur
    L_deltaT = []     # Différences de température par segment
    L_Longueur = []   # Longueurs associées

    compteur = 1      # Épaisseur cumulée du matériau traversé
    T_0 = T0[0, i]    # Température côté extérieur
    T_lim = T0[0, i]  # Limite initiale pour les différences de température

    # --- Calcul des flux linéiques du mur de référence ---
    while i < geo.n - 1:  # Jusqu'à l'air intérieur chauffé
        # Si changement de matériau, calcul du flux pour le segment précédent
        if geo0.G[1, i + 1, 0] != geo0.G[1, i, 0]:
            L0i = compteur * geo.dx

            # Ajout du flux linéique (analogie électrique)
            L_U0j.append(geo0.G[1, i, 1] / L0i)
            L_deltaT.append(T0[0, i] - T_lim)
            L_Longueur.append(L0i)

            # Réinitialisation pour le nouveau matériau
            compteur = 1
            T_lim = T0[0, i]
        else:
            compteur += 1
        i += 1

    # Conversion en arrays numpy pour faciliter les calculs vectorisés
    L_U0j = np.array(L_U0j)
    L_deltaT = np.array(L_deltaT)
    L_Longueur = np.array(L_Longueur)

    # --- Détection des indices caractéristiques ---
    if geo.Nom_Geometrie == 'plancher simple':
        p = int((geo.li + geo.e) / geo.dx + 1)
        q1 = int((geo.hi - geo.eps_iso) / geo.dy)
        q2 = int((geo.hi + geo.h) / geo.dy + 1)
    else:
        p = int((geo.li + geo.e) / geo.dx + 1)
        q1 = int((geo.hi - geo.eps_iso) / geo.dy)
        q2 = int((geo.hi + geo.h + geo.eps_iso) / geo.dy + 1)

    T_lim = T0[0, i]  # Température côté intérieur

    # --- Calcul des flux linéiques ---
    Phi = geo.P * (
        np.sum(J[:q1, p, 0]) * geo.dy / geo.hi +
        np.sum(J[q2:, p, 0]) * geo.dy / geo.hi -
        np.sum(np.abs(J[q1 + 1, p:, 1])) * geo.dx / geo.Lmax -
        np.sum(np.abs(J[q2 - 1, p:, 1])) * geo.dx / geo.Lmax
    )

    # Flux de référence (mur seul)
    Phi_default = -np.sum(L_U0j * L_deltaT) * geo.P

    # --- Calcul de Psi ---
    Psi = -(Phi - Phi_default) / (T_lim - T_0)

    return Psi, p, q1, q2

