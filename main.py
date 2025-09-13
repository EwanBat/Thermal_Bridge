### Import des librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from geometrie import Geometry
from fonction_calcul import schema_jacobi  # si dans le même package
from fonction_therm import CL, Jth, Psi_lin
from optimisation import optimisation_floor, optimisation_planelle

plt.close('all')

## Simulation d'un pont thermique en 2D
# Paramétrage
################################################################
# PARAMÈTRES DE SIMULATION
################################################################

# Paramètres économiques
prix = 0.2062e-3  # Prix du kWh en euros

# Paramètres géométriques du bâtiment
e = 0.22      # Épaisseur du mur en béton (m)
h = 0.2       # Épaisseur du plancher (m) (optimal: 0.12m)
li = 0.1      # Épaisseur de l'isolant (m) (optimal: 0.07m)
eps_iso = 0   # Épaisseur du plancher isolant (m)
Coef = 0.6   # Position relative planelle/rupteur (optimal: 8/10)

# Dimensions du domaine de calcul
Lmax = 0.8    # Largeur du domaine (m)
Hmax = 0.8    # Hauteur du domaine (m)
P = 2         # Profondeur de la pièce (m)
S = P*Hmax    # Surface de la paroi (m²)
hi = round((Hmax-h)/2, 3)  # Demi-hauteur de l'isolant (m)

# Propriétés thermiques des matériaux
## Béton
Cth1 = 1.65   # Conductivité thermique (W/m.K)
Cp1 = 1000    # Capacité thermique (J/kg.K)
rho1 = 2150   # Masse volumique (kg/m³)

## Isolant
Cth2 = 0.03   # Conductivité thermique (W/m.K)

## Air
Cth0 = 0.024  # Conductivité thermique (W/m.K)
h0 = 6        # Coefficient d'échange surfacique intérieur (W/m².K)
hm1 = 15      # Coefficient d'échange surfacique extérieur (W/m².K)

## Acier (pour rupteur)
Cth3 = 50.2   # Conductivité thermique (W/m.K)

# Paramètres thermiques
Tint = 19 + 273    # Température intérieure (K)
Text1 = 10 + 273   # Température extérieure initiale (K)
Text2 = -10 + 273  # Température extérieure finale (K)
a = Cth1/Cp1/rho1  # Diffusivité thermique du béton (m²/s)

# Paramètres numériques
n, m = 100, 100    # Nombre de points de discrétisation en x et y
dx = Lmax/n        # Pas d'espace en x (m)
dy = Hmax/m        # Pas d'espace en y (m)
dt = 0.25*min(dx,dy)**2/a  # Pas de temps (s) - Condition CFL
tmax = 1500        # Durée totale de simulation (s)
w = 2/(1 + np.pi/n)  # Coefficient de relaxation

## Fonctions utilitaires
x,y = np.arange(0,Lmax,dx),np.arange(0,Hmax,dy)
L_N,L_M = [i for i in range(n)],[j for j in range(m)]
ytick_label = ['Extérieur','Air intérieur','Béton','Matériau isolant']

def discrete_matshow(data, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    Cmap = plt.get_cmap('rainbow', 4) # type: ignore
    mat = ax.matshow(data, cmap=Cmap, vmin=-1 - 0.5, vmax=2 + 0.5)
    cax = plt.colorbar(mat, ticks=np.arange(-1, 2+1), ax=ax)
    cax.ax.set_yticklabels(ytick_label)
    ax.set_title('Géométrie ' + str(geo.Nom_Geometry))
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
    ax2.set_title('T(x,y) pour ' + geo.Nom_Geometry)
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

    print('Le coefficient linéaire de pont thermique est ', Psi, ' W/(m.K) (Géométrie ' + geo.Nom_Geometry + ') \n')
    print('Le coût de la perte due au pont thermique est', prix*24*(Tint-Text2)*P*Psi, '€ en 24 heures \n')
    plt.show()

###################################### Commande à l'utilisateur #########################
print("Choose which activity you want to do : \n")
print("1 : Simple simulation of a geometry \n")
print("2 : Optimization of the insulation thickness on the floor \n")
print("3 : Optimization of the position of the thermal bridge breaker or the insulation strip \n")
answer = input("Your choice (1, 2 or 3) : ")

geo0 = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, P, Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef) # Geométrie de référence : mur
geo0.mur()
################################################################### Application #####################################################################
if answer == '1':
    geo = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, P, Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef)
    geo.rupteur()

    # T,Temps_stat = Initialisation(),0
    T0, Temps_stat = schema_jacobi(geo0,CL)
    T, Temps_stat = schema_jacobi(geo,CL)

    print('Temps prit pour aller jusqu au régime stationnaire :',Temps_stat,' min \n')

    J = Jth(geo,T)

    Psi, p, q1, q2 = Psi_lin(geo0,T0,geo,T,J)
    Psi = round(Psi,4)

    Affiche(geo,T0,T,J,Psi)

###################################### Réduction de Psi au maximum ###########################################
################################ Optimisation de l'épaisseur de l'isolant ################################
elif answer == '2':
    ''' Evolution de l'épaisseur d'isolant '''

    geo = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, P, Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef)
    method_name = 'planelle_isolant'

    li_min, li_max, n_steps = 0.05, 0.22, 5
    best_li, best_psi, L_li, L_psi = optimisation_floor(geo0, geo, li_min, li_max, n_steps, method_name)
    print(f"Épaisseur optimale de l'isolant : {best_li:.4f} m")
    print(f"Coefficient linéique Ψ correspondant : {best_psi:.4f} W/m·K")

    plt.figure('Optimisation de l\'épaisseur de l\'isolant')
    plt.plot(L_li, L_psi, marker='o')
    plt.title('Optimisation de l\'épaisseur de l\'isolant')
    plt.xlabel('Épaisseur de l\'isolant (m)')
    plt.ylabel('Coefficient linéique Ψ (W/m·K)')
    plt.grid()
    plt.show()

elif answer == '3':
    ''' Emplacement de la planelle ou du rupteur '''

    geo = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, P, Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef)
    method_name = "rupteur"  # ou 'planelle_isolant'
    best_Coef, best_psi, L_Coef, L_psi = optimisation_planelle(geo0, geo, method_name)
    print(f"Position optimale de la planelle/rupteur : {best_Coef:.2f} m")
    print(f"Coefficient linéique Ψ correspondant : {best_psi:.4f} W/(m·K)")

    plt.figure('Optimisation de la position de la planelle/rupteur')
    plt.plot(L_Coef, L_psi, marker='o')
    plt.title('Optimisation de la position de la planelle/rupteur')
    plt.xlabel('Position de la planelle/rupteur (m)')
    plt.ylabel('Coefficient linéique Ψ (W/m·K)')
    plt.grid()
    plt.show()


# for g in L_Geo[2:4]:
#     L_Psi = []
#     L_Temps = []
#     print(str(g)+'\n')
#     for i in L_Coef:
#         Coef = i/10
#         print(Coef)

#         G0 = Geometry_mur(n,m)
#         G = g(n,m)

#         T0,Temps_stat = schema_jacobi(G0,CL)
#         T,Temps_stat = schema_jacobi(G,CL)

#         J = Jth(G,T)

#         Psi = round(Psi_lin(G0,T0,G,T,J),4)
#         L_Psi.append(Psi)
#         L_Temps.append(Temps_stat)
#     Moy_g = Moyenne(L_Psi)
#     Dico_Moyenne['G'+Nom_Geometry] = Moy_g

#     plt.figure('Psi selon planelle/rupteur')
#     plt.plot(L_Coef*e/10,L_Psi,label=Nom_Geometry)

#     plt.figure('Temps selon planelle/rupteur')
#     plt.plot(L_Coef*e/10,L_Temps,label=Nom_Geometry)

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
