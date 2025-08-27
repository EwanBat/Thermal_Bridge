## Import des librairies
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from Geometry import Geometrie
import matplotlib.gridspec as gridspec

plt.close('all')

## Simulation d'un pont thermique en 2D
# Paramétrage
prix = 0.2062e-3 #Prix du Watt.heure en €

e = 0.22 #Epaisseur du béton en m sur x
Coef = 5/10 #Coefficient pour l'emplacement de la planelle ou du rupteur (5/10 par défaut ; 8/10 optimal)
h = 0.2 #Epaisseur du sol en m sur y (0.2 par défaut ; 0.12 optimal)

li = 0.1 #Epaisseur de l'isolant en m sur x (0.1 par défaut ; 0.07 optimal)
ep_iso = 0 #Epaisseur du plancher isolant en m sur x

Lmax,Hmax = 0.8,0.8 #Dimension de l'environnement en m
tmax = 1500 #Durée de l'expérience en s à partir de t=0
hi = round((Hmax-h)/2,3) #Demi hauteur de l'isolant en m
P = 2 #Profondeur de la pièce en m
S = P*Hmax

Cth1 = 1.65 #Conductivité thermique du mur en W/K/m
Cp1 = 1000 #Capacité thermique du mur en J/kg/K
rho1 = 2150 #Masse volumique du mur en kg/m3

Cth2 = 0.03 #Conductivité thermique de l'isolant en W/K/m

h0 = 6 #Coefficient de transfert thermique de l'air intérieur en W/m2/K
Cth0 = 0.024 #Conductivité thermique de l'air intérieur en W/K/m

hm1 = 15 #Coefficient de transfert thermique de l'air extérieur en W/m2/K

Cth3 = 50.2 #Conductivité thermique de l'acier en W/K/m

a = Cth1/Cp1/rho1

Tint = 19 + 273 #Température de la pièce en K
Text1 = 10 + 273 #Température extérieure pour t<0 en K
Text2 = -10 + 273 #Température extérieure pour t>0 en K

n,m = 100,100 #Nombres de points de discrétisation sur x et y
dx,dy = Lmax/n,Hmax/m
dt = 0.25*min(dx,dy)**2/a
w = 2/(1 + np.pi/n)
x,y = np.arange(0,Lmax,dx),np.arange(0,Hmax,dy)
L_N,L_M = [i for i in range(n)],[j for j in range(m)]

ytick_label = ['Extérieur','Air intérieur','Béton','Matériau isolant']
## Fonctions utilitaires
def discrete_matshow(data, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    Cmap = plt.get_cmap('rainbow', 4)
    mat = ax.matshow(data, cmap=Cmap, vmin=-1 - 0.5, vmax=2 + 0.5)
    cax = plt.colorbar(mat, ticks=np.arange(-1, 2+1), ax=ax)
    cax.ax.set_yticklabels(ytick_label)
    ax.set_title('Géométrie ' + str(geo.Nom_Geometrie))
    ax.set_xlabel('x-axis (m)')
    ax.set_ylabel('y-axis (m)')
    ax.set_xticks(L_N[0::n//10]+L_N[-1:])
    ax.set_xticklabels(list(x[0::n//10])+[Lmax])
    ax.set_yticks(L_M[0::n//10]+L_M[-1:])
    ax.set_yticklabels(list(y[0::n//10])+[Hmax])

def Affiche(geo, T0, T, J, Psi):
    fig = plt.figure(constrained_layout=True, figsize=(15, 5))
    gs = gridspec.GridSpec(1, 3, figure=fig)

    # Géométrie
    ax0 = fig.add_subplot(gs[0, 0])
    discrete_matshow(np.round(geo.G[:,:,0], 2), ax=ax0)

    # Température mur
    ax1 = fig.add_subplot(gs[0, 1])
    c1 = ax1.contourf(np.round(T0, 2), 200, cmap='jet')
    fig.colorbar(c1, ax=ax1, label='Température (°C)')
    ax1.set_title('T(x,y) pour un mur')
    ax1.set_xlabel('x-axis (m)')
    ax1.set_ylabel('y-axis (m)')
    ax1.set_xticks(L_N[0::n//10]+L_N[-1:])
    ax1.set_xticklabels(list(x[0::n//10])+[Lmax])
    ax1.set_yticks(L_M[0::n//10]+L_M[-1:])
    ax1.set_yticklabels(list(y[0::n//10])+[Hmax])

    # Température géométrie
    ax2 = fig.add_subplot(gs[0, 2])
    c2 = ax2.contourf(np.round(T, 2), 200, cmap='jet')
    fig.colorbar(c2, ax=ax2, label='Température (°C)')
    ax2.set_title('T(x,y) pour ' + geo.Nom_Geometrie)
    ax2.set_xlabel('x-axis (m)')
    ax2.set_ylabel('y-axis (m)')
    ax2.set_xticks(L_N[0::n//10]+L_N[-1:])
    ax2.set_xticklabels(list(x[0::n//10])+[Lmax])
    ax2.set_yticks(L_M[0::n//10]+L_M[-1:])
    ax2.set_yticklabels(list(y[0::n//10])+[Hmax])

    # Ajout des vecteurs de flux thermique sur la dernière figure
    Compteur = 0
    for i in range(n):
        for j in range(m):
            if i % 10 == 0 and j % 10 == 0:
                Grad = J[j, i]
                di, dj = Grad[0], Grad[1]
                if Compteur == 0:
                    ax2.arrow(i, j, di, dj, label='Jth = g(x,t)', width=0.3, color='black')
                    Compteur = 1
                else:
                    ax2.arrow(i, j, di, dj, width=0.3, color='black')

    # Affichage des lignes de séparation (si p, q1, q2 existent)
    try:
        ax2.plot([p, n], [q1, q1], color='white')
        ax2.plot([p, n], [q2, q2], color='white')
        ax2.plot([p, p], [0, q1], color='white')
        ax2.plot([p, p], [q2, m], color='white')
    except NameError:
        pass

    ax2.legend()

    print('Le coefficient linéaire de pont thermique est ', Psi, ' W/(m.K) (Géométrie ' + geo.Nom_Geometrie + ') \n')
    print('Le coût de la perte due au pont thermique est', prix*24*(Tint-Text2)*P*Psi, '€ en 24 heures \n')
    plt.show()

def gradxy(T,i,j):
    if i == n-1:
        if j == m-1:
            Y,X = (T[m-2,i]-T[j,i])/dy,(T[j,n-2]-T[j,n-1])/dx
        else:
            Y,X = (T[j+1,i]-T[j,i])/dy,(T[j,n-2]-T[j,n-1])/dx
    else:
        if j == m-1:
            Y,X = (T[m-2,i]-T[j,i])/dy,(T[j,i+1]-T[j,i])/dx
        else:
            Y,X = (T[j+1,i]-T[j,i])/dy,(T[j,i+1]-T[j,i])/dx
    return np.array([X,Y])

def Moyenne(L):
    Somme = 0
    for l in L:
        Somme += l
    return Somme/len(L)

def Initialisation(): #Initialise de manière affine la température en K
    T = np.zeros([m,n])
    for j in range(m):
        T[j,:] = [Text1 - (Text1-Tint)*dx*i/Lmax for i in range(n)]
    return T

def calc(T,i,j): #Calcul du terme en j,i par méthode de différence finie
    Terme = w*(1/4)*(T[j+1,i] + T[j-1,i] + T[j,i+1] + T[j,i-1]) + (1-w)*T[j,i]
    return Terme

def itere(geo, T, Condition): #Itère le calcul par différence finie selon si j,i est en condition limite ou non
    Somme = 0
    for i in range(n):
        for j in range(m):
            Bool = Condition(geo, T, i, j) #Sert à savoir si on est en condition limite par le biai d'un booléen
            if not Bool[0]: #Si on est pas en condition limite on fait le calcul par la méthode de différence finie
                Terme = calc(T, i, j)
                Somme += (Terme - T[j,i])**2
                T[j,i] = Terme
            else: #Sinon le calcul est déjà fait par la fonction CL donc on l'utilise
                Terme = Bool[1]
                Somme += (Terme - T[j,i])**2
                T[j,i] = Terme
    return Somme #Renvoie une somme pour savoir l'écart entre l'état de départ et le régime permanent

def schema_jacobi(geo, Condition):
    T = Initialisation()
    Ind = itere(geo,T,Condition)
    Compteur = 1
    while Ind >= 1e-1: #On recommence jusqu'à arriver en régime permanent
        Ind = itere(geo,T,Condition)
        Compteur += 1
    Temps_stat = round(Compteur*dt/60,1) # Temps pour arriver au régime stationnaire en min
    T = T - 273
    return T,Temps_stat

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
        Tij = Text2
        return True, Tij

    # --- BORD DROIT : Conditions symétriques ou température intérieure ---
    elif i == n - 1:
        if j == 0:  # Coin supérieur droit
            if geo.G[j, i][0] != 0:  # Non-air → symétrie
                Tij = w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j + 1, i] + T[m - 2, i]) \
                      + (1 - w) * T[j, i]
            else:  # Air → température intérieure
                Tij = Tint

        elif j == m - 1:  # Coin inférieur droit
            if geo.G[j, i][0] != 0:  # Non-air → symétrie
                Tij = w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j - 1, i] + T[1, i]) \
                      + (1 - w) * T[j, i]
            else:  # Air → température intérieure
                Tij = Tint

        else:  # Bord droit sans être dans un coin
            if geo.G[j, i][0] != 0:
                Tij = w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j + 1, i] + T[j - 1, i]) \
                      + (1 - w) * T[j, i]
            else:
                Tij = Tint

        return True, Tij

    # --- BORD SUPÉRIEUR ---
    elif j == 0:
        Tij = w * (1 / 4) * (T[j, i + 1] + T[j, i - 1] + T[j + 1, i] + T[m - 2, i]) \
              + (1 - w) * T[j, i]
        return True, Tij

    # --- BORD INFÉRIEUR ---
    elif j == m - 1:
        Tij = w * (1 / 4) * (T[j, i + 1] + T[j, i - 1] + T[j - 1, i] + T[1, i]) \
              + (1 - w) * T[j, i]
        return True, Tij

    # --- INTERFACE ENTRE DEUX MATÉRIAUX SUR L'AXE X ---
    elif geo.G[j, i][0] != geo.G[j, i + 1][0] or geo.G[j, i][0] != geo.G[j, i - 1][0]:
        lamb1 = geo.G[j, i - 1][1]
        lamb2 = geo.G[j, i + 1][1]

        # Cas avec convection (type -1)
        if geo.G[j, i - 1][0] == -1 and geo.G[j, i][0] != geo.G[j, i - 1][0]:
            Tij = (lamb2 * T[j, i + 1] / dx + T[j, i - 1] * hm1) / (lamb2 / dx + hm1)
        elif geo.G[j, i + 1][0] == -1 and geo.G[j, i][0] != geo.G[j, i + 1][0]:
            Tij = (lamb2 * T[j, i - 1] / dx + T[j, i + 1] * hm1) / (lamb2 / dx + hm1)
        else:  # Interface purement conductrice
            Tij = (lamb1 * T[j, i - 1] + lamb2 * T[j, i + 1]) / (lamb2 + lamb1)

        return True, Tij

    # --- INTERFACE ENTRE DEUX MATÉRIAUX SUR L'AXE Y ---
    elif geo.G[j, i][0] != geo.G[j + 1, i][0] or geo.G[j, i][0] != geo.G[j - 1, i][0]:
        lamb1 = geo.G[j - 1, i][1]
        lamb2 = geo.G[j + 1, i][1]

        # Cas avec convection (type -1)
        if geo.G[j - 1, i][0] == -1 and geo.G[j, i][0] != geo.G[j - 1, i][0]:
            Tij = (lamb2 * T[j + 1, i] / dx + T[j - 1, i] * hm1) / (lamb2 / dx + hm1)
        elif geo.G[j + 1, i][0] == -1 and geo.G[j, i][0] != geo.G[j + 1, i][0]:
            Tij = (lamb2 * T[j - 1, i] / dx + T[j + 1, i] * hm1) / (lamb2 / dx + hm1)
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
    J = np.zeros((m, n, 2))

    # Parcours de chaque point de la grille
    for i in range(n):       # axe x
        for j in range(m):   # axe y
            # Gradient local de température (∇T)
            Grad = gradxy(T, i, j)

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
    while i < n - 1:  # Jusqu'à l'air intérieur chauffé
        # Si changement de matériau, calcul du flux pour le segment précédent
        if geo0.G[1, i + 1, 0] != geo0.G[1, i, 0]:
            L0i = compteur * dx

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
    global p, q1, q2  # Indices des interfaces mur/isolation et air intérieur
    if geo.Nom_Geometrie == 'plancher simple':
        p = int((li + e) / dx + 1)
        q1 = int((hi - ep_iso) / dy)
        q2 = int((hi + h) / dy + 1)
    else:
        p = int((li + e) / dx + 1)
        q1 = int((hi - ep_iso) / dy)
        q2 = int((hi + h + ep_iso) / dy + 1)

    T_lim = T0[0, i]  # Température côté intérieur

    # --- Calcul des flux linéiques ---
    Phi = P * (
        np.sum(J[:q1, p, 0]) * dy / hi +
        np.sum(J[q2:, p, 0]) * dy / hi -
        np.sum(np.abs(J[q1 + 1, p:, 1])) * dx / Lmax -
        np.sum(np.abs(J[q2 - 1, p:, 1])) * dx / Lmax
    )

    # Flux de référence (mur seul)
    Phi_default = -np.sum(L_U0j * L_deltaT) * P

    print(Phi, Phi_default, '\n')

    # --- Calcul de Psi ---
    Psi = -(Phi - Phi_default) / (T_lim - T_0)

    return Psi


################################################################### Application #####################################################################

geo0 = Geometrie(n,m,e,li,hi,h,dx,dy,Cth0,Cth1,Cth2)
geo0.mur()
geo = Geometrie(n,m,e,li,hi,h,dx,dy,Cth0,Cth1,Cth2)
geo.planelle_isolant()

# T,Temps_stat = Initialisation(),0
T0,Temps_stat = schema_jacobi(geo0,CL)
T,Temps_stat = schema_jacobi(geo,CL)

print('Temps prit pour aller jusqu au régime stationnaire :',Temps_stat,' min \n')

J = Jth(geo,T)

Psi = round(Psi_lin(geo0,T0,geo,T,J),4)

Affiche(geo,T0,T,J,Psi)

############################################################# Réduction de Psi au maximum ###########################################################################
''' Evolution de l'épaisseur d'isolant '''
L_li = np.round(np.arange(0.05,0.21,0.01),4)
L_h = np.arange(0.12,0.21,0.01)
L_Coef = np.array([i for i in range(1,9)])
""" L_Geo = [Geometrie_ss,Geometrie_classique,Geometrie_rupteur,Geometrie_planelle_isolant,Geometrie_classique_plancher_isolant]
L_Geo = L_Geo
Dico_Moyenne = {} """

# for g in L_Geo:
#     L_Psi = []
#     L_Temps = []
#     print(str(g)+'\n')
#     for li in list(L_li):
#         li = round(li,3)
#         print(li)

#         G0 = Geometrie_mur(n,m)
#         G = g(n,m)

#         T0,Temps_stat = schema_jacobi(G0,CL)
#         T,Temps_stat = schema_jacobi(G,CL)

#         J = Jth(G,T)

#         Psi = round(Psi_lin(G0,T0,G,T,J),4)
#         L_Psi.append(Psi)
#         L_Temps.append(Temps_stat)

#     Moy_g = Moyenne(L_Psi)
#     Dico_Moyenne['G '+Nom_Geometrie] = Moy_g

#     plt.figure('Psi selon isolant')
#     plt.plot(L_li,L_Psi,label=Nom_Geometrie)

#     plt.figure('Temps selon isolant')

#     plt.plot(L_li,L_Temps,label=Nom_Geometrie)

#     print('\n')

# plt.figure('Efficacité selon isolant')
# plt.title('Psi moyen selon isolant')
# plt.xylabel('Psi (W/(m.k))')
# plt.bar([i for i in range(len(L_Geo))],Dico_Moyenne.values())
# plt.xticks([i for i in range(len(L_Geo))],Dico_Moyenne.keys())

# plt.figure('Psi selon isolant')
# plt.title('Psi = h(li), li : épaisseur du matériau isolant')
# plt.xlabel('épaisseur du matériau isolant (m)')
# plt.ylabel('Psi (W/(m.K))')
# plt.legend()

# plt.figure('Temps selon isolant')
# plt.title('Time = k(li), li : épaisseur du matériau isolant')
# plt.xlabel('épaisseur du matériau isolant (m)')
# plt.ylabel('Temps mis pour atteindre l état stationnaire (min)')
# plt.legend()

# plt.show()

''' Evolution de l'épaisseur de la dalle de béton '''

# for g in L_Geo:
#     L_Psi = []
#     L_Temps = []
#     print(str(g)+'\n')
#     for h in L_h:
#         h = round(h,4)
#         hi = round((Hmax-h)/2,3) #Demi hauteur de l'isolant en m
#
#         print(h)
#
#         G0 = Geometrie_mur(n,m)
#         G = g(n,m)
#
#         T0,Temps_stat = schema_jacobi(G0,CL)
#         T,Temps_stat = schema_jacobi(G,CL)
#
#         J = Jth(G,T)
#
#         Psi = round(Psi_lin(G0,T0,G,T,J),4)
#         L_Psi.append(Psi)
#         L_Temps.append(Temps_stat)
#
#     Moy_g = Moyenne(L_Psi)
#
#     Dico_Moyenne['G'+Nom_Geometrie] = Moy_g
#     plt.figure('Psi selon chape')
#     plt.plot(L_h,L_Psi,label=Nom_Geometrie)
#
#     plt.figure('Temps selon chape')
#     plt.plot(L_h,L_Temps,label=Nom_Geometrie)
#
#     print('\n')
#
# plt.figure('Efficacité selon chape')
# plt.title('Psi moyen selon chape')
# plt.xlabel('Psi (W/(m.k))')
# plt.ylabel('épaisseur de la chape (m)')
# plt.bar([i for i in range(len(L_Geo))],Dico_Moyenne.values())
# plt.xticks([i for i in range(len(L_Geo))],Dico_Moyenne.keys())
#
# plt.figure('Psi selon chape')
# plt.title('Psi = h(h), h : épaisseur de la chape de béton')
# plt.xlabel('épaisseur de la chape (m)')
# plt.ylabel('Psi (W/(m.K))')
# plt.legend()
#
# plt.figure('Temps selon chape')
# plt.title('Time = k(h), h : épaisseur de la chape de béton')
# plt.xlabel('épaisseur de la chape (m)')
# plt.ylabel('Temps mis pour atteindre l état stationnaire (min)')
# plt.legend()
#
# plt.show()

''' Emplacement de la planelle ou du rupteur '''

# for g in L_Geo[2:4]:
#     L_Psi = []
#     L_Temps = []
#     print(str(g)+'\n')
#     for i in L_Coef:
#         Coef = i/10
#         print(Coef)

#         G0 = Geometrie_mur(n,m)
#         G = g(n,m)

#         T0,Temps_stat = schema_jacobi(G0,CL)
#         T,Temps_stat = schema_jacobi(G,CL)

#         J = Jth(G,T)

#         Psi = round(Psi_lin(G0,T0,G,T,J),4)
#         L_Psi.append(Psi)
#         L_Temps.append(Temps_stat)
#     Moy_g = Moyenne(L_Psi)
#     Dico_Moyenne['G'+Nom_Geometrie] = Moy_g

#     plt.figure('Psi selon planelle/rupteur')
#     plt.plot(L_Coef*e/10,L_Psi,label=Nom_Geometrie)

#     plt.figure('Temps selon planelle/rupteur')
#     plt.plot(L_Coef*e/10,L_Temps,label=Nom_Geometrie)

#     print('\n')

# plt.figure('Efficacité selon planelle / rupteur')
# plt.title('Psi moyen selon emplacement')
# plt.xlabel('Psi (W/(m.k))')
# plt.ylabel('emplacement de la planelle / du rupteur (m)')
# plt.bar([i for i in range(len(L_Geo[2:4]))],Dico_Moyenne.values())
# plt.xticks([i for i in range(len(L_Geo[2:4]))],Dico_Moyenne.keys())

# plt.figure('Psi selon planelle/rupteur')
# plt.title('Psi = h(Coef), Coef : Emplacement de la planelle/rupteur')
# plt.xlabel('Emplacement de la planelle/rupteur (m)')
# plt.ylabel('Psi (W/(m.K))')
# plt.legend()

# plt.figure('Temps selon planelle/rupteur')
# plt.title('Time = k(h), h : Emplacement de la planelle/rupteur')
# plt.xlabel('Emplacement de la planelle/rupteur (m)')
# plt.ylabel('Temps mis pour atteindre l état stationnaire (min)')
# plt.legend()

# plt.show()
