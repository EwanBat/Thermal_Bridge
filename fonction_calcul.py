import numpy as np

def gradxy(geo, T, i, j):
    """
    Calcule le gradient de température aux coordonnées (i,j)
    
    Args:
        geo (Geometrie): Objet contenant les paramètres géométriques
        T (np.array): Matrice des températures
        i,j (int): Indices de position dans la grille
        
    Returns:
        np.array: Gradient [dT/dx, dT/dy]
    """
    if i == geo.n-1:  # Bord droit
        if j == geo.m-1:  # Coin inférieur droit
            Y = (T[geo.m-2,i]-T[j,i])/geo.dy  # Différence finie arrière en y
            X = (T[j,geo.n-2]-T[j,geo.n-1])/geo.dx  # Différence finie arrière en x
        else:  # Bord droit (non coin)
            Y = (T[j+1,i]-T[j,i])/geo.dy  # Différence finie avant en y
            X = (T[j,geo.n-2]-T[j,geo.n-1])/geo.dx  # Différence finie arrière en x
    else:  # Points intérieurs ou bord gauche
        if j == geo.m-1:  # Bord inférieur
            Y = (T[geo.m-2,i]-T[j,i])/geo.dy  # Différence finie arrière en y
            X = (T[j,i+1]-T[j,i])/geo.dx  # Différence finie avant en x
        else:  # Points intérieurs
            Y = (T[j+1,i]-T[j,i])/geo.dy  # Différence finie avant en y
            X = (T[j,i+1]-T[j,i])/geo.dx  # Différence finie avant en x
    return np.array([X,Y])

def Initialisation(geo):
    """
    Initialise le champ de température avec une variation linéaire
    entre Text1 (extérieur) et Tint (intérieur)
    
    Args:
        geo (Geometrie): Objet contenant les paramètres géométriques
        
    Returns:
        np.array: Matrice des températures initiales
    """
    T = np.zeros([geo.m, geo.n])
    for j in range(geo.m):
        T[j,:] = [geo.Text1 - (geo.Text1-geo.Tint)*geo.dx*i/geo.Lmax for i in range(geo.n)]
    return T

def calc(geo, T, i, j):
    """
    Calcule la nouvelle température au point (i,j) par la méthode de relaxation
    
    Args:
        geo (Geometrie): Objet contenant les paramètres géométriques
        T (np.array): Matrice des températures actuelles
        i,j (int): Indices de position dans la grille
        
    Returns:
        float: Nouvelle température calculée
    """
    # Méthode de relaxation avec coefficient w
    Terme = geo.w*(1/4)*(T[j+1,i] + T[j-1,i] + T[j,i+1] + T[j,i-1]) + (1-geo.w)*T[j,i]
    return Terme

def itere(geo, T, Condition):
    """
    Effectue une itération complète sur toute la grille
    
    Args:
        geo (Geometrie): Objet contenant les paramètres géométriques
        T (np.array): Matrice des températures
        Condition (function): Fonction définissant les conditions aux limites
        
    Returns:
        float: Somme des écarts quadratiques entre deux itérations
    """
    Somme = 0
    for i in range(geo.n):
        for j in range(geo.m):
            Bool, val = Condition(geo, T, i, j)  # Vérifie si point aux conditions limites
            if not Bool:  # Point interne : calcul par relaxation
                Terme = calc(geo, T, i, j)
            else:  # Point limite : utilise la valeur imposée
                Terme = val
            Somme += (Terme - T[j,i])**2  # Accumule l'écart
            T[j,i] = Terme  # Met à jour la température
    return Somme

def schema_jacobi(geo, Condition):
    """
    Résout l'équation de la chaleur par la méthode de Jacobi jusqu'à convergence
    
    Args:
        geo (Geometrie): Objet contenant les paramètres géométriques
        Condition (function): Fonction définissant les conditions aux limites
        
    Returns:
        tuple: (Matrice des températures finale, Temps de simulation)
    """
    T = Initialisation(geo)  # Initialisation du champ
    Ind = itere(geo,T,Condition)  # Première itération
    Compteur = 1
    
    # Itère jusqu'à convergence (seuil 0.1)
    while Ind >= 1e-1:
        Ind = itere(geo,T,Condition)
        Compteur += 1
        
    # Calcul du temps de simulation et conversion en Celsius
    Temps_stat = round(Compteur*geo.dt/60,1)
    T = T - 273
    
    return T,Temps_stat