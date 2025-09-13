### Library imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from geometrie import Geometry
from fonction_calcul import schema_jacobi
from fonction_therm import CL, Jth, Psi_lin
from optimisation import optimisation_floor, optimisation_planelle

plt.close('all')

## 2D thermal bridge simulation
# Settings
################################################################
# SIMULATION PARAMETERS
################################################################

# Economic parameters
prix = 0.2062e-3  # Price per kWh in euros

# Building geometric parameters
e = 0.22      # Concrete wall thickness (m)
h = 0.2       # Floor thickness (m) (optimal: 0.12m)
li = 0.1      # Insulation thickness (m) (optimal: 0.07m)
eps_iso = 0.1 # Floor insulation thickness (m)
Coef = 0.6    # Relative position of block/thermal break (optimal: 8/10)

# Calculation domain dimensions
Lmax = 0.8    # Domain width (m)
Hmax = 0.8    # Domain height (m)
P = 2         # Room depth (m)
S = P*Hmax    # Wall surface area (m²)
hi = round((Hmax-h)/2, 3)  # Height of interface air-wall (m)

# Thermal properties of materials
## Concrete
Cth1 = 1.65   # Thermal conductivity (W/m.K)
Cp1 = 1000    # Heat capacity (J/kg.K)
rho1 = 2150   # Density (kg/m³)

## Insulation
Cth2 = 0.03   # Thermal conductivity (W/m.K)

## Air
Cth0 = 0.024  # Thermal conductivity (W/m.K)
h0 = 6        # Interior surface heat transfer coefficient (W/m².K)
hm1 = 15      # Exterior surface heat transfer coefficient (W/m².K)

## Steel (for thermal break)
Cth3 = 50.2   # Thermal conductivity (W/m.K)

# Thermal parameters
Tint = 19 + 273    # Interior temperature (K)
Text1 = 10 + 273   # Initial exterior temperature (K)
Text2 = -10 + 273  # Final exterior temperature (K)
a = Cth1/Cp1/rho1  # Thermal diffusivity of concrete (m²/s)

# Numerical parameters
n, m = 100, 100    # Number of discretization points in x and y
dx = Lmax/n        # Space step in x (m)
dy = Hmax/m        # Space step in y (m)
dt = 0.25*min(dx,dy)**2/a  # Time step (s) - CFL condition
tmax = 1500        # Total simulation duration (s)
w = 2/(1 + np.pi/n)  # Relaxation coefficient

## Utility functions
x,y = np.arange(0,Lmax,dx),np.arange(0,Hmax,dy)
L_N,L_M = [i for i in range(n)],[j for j in range(m)]
ytick_label = ['Exterior','Interior air','Concrete','Insulation material']

def discrete_matshow(data, ax=None):
    """
    Display geometry with material colors and labels
    
    Parameters:
        data: Material type matrix from geo.G[:,:,0]
        ax: Matplotlib axis for plotting
    """
    if ax is None:
        fig, ax = plt.subplots()
        
    # Define material colors and labels
    materials = {
        -1: ('lightblue', 'External air'),
        0: ('white', 'Interior air'),
        1: ('gray', 'Concrete'),
        2: ('yellow', 'Insulation')
    }
    
    # Create custom colormap
    colors = [materials[i][0] for i in [-1, 0, 1, 2]]
    cmap = plt.matplotlib.colors.ListedColormap(colors)  # type: ignore
    
    # Display material matrix
    mat = ax.matshow(np.flip(data, axis=0), cmap=cmap, vmin=-1.5, vmax=2.5)
    
    # Create custom legend
    patches = [plt.matplotlib.patches.Patch(color=materials[i][0],  # type: ignore
                                          label=materials[i][1]) 
              for i in [-1, 0, 1, 2]]
    ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Set labels and title
    ax.set_title('Geometry: ' + str(geo.Nom_Geometry))
    ax.set_xlabel('x-axis (m)')
    ax.set_ylabel('y-axis (m)')
    
    # Set ticks
    ax.set_xticks(L_N[0::n//10]+L_N[-1:])
    ax.set_xticklabels([f'{x:.2f}' for x in np.arange(0, Lmax+dx, Lmax/10)])
    ax.set_yticks(L_M[0::n//10]+L_M[-1:])
    ax.set_yticklabels([f'{y:.2f}' for y in np.arange(0, Hmax+dy, Hmax/10)])

def Affiche(geo, T0, T, J, Psi, p, q1, q2):
    """
    Display geometry, temperature fields and heat flux with boundary lines
    
    Parameters:
        geo: Geometry object
        T0: Reference temperature field
        T: Temperature field with thermal bridge
        J: Heat flux field
        Psi: Linear thermal transmittance
        p, q1, q2: Boundary indices for wall and floor
    """
    fig = plt.figure(constrained_layout=True, figsize=(15, 5))
    gs = gridspec.GridSpec(1, 3, figure=fig)

    # Geometry with materials
    ax0 = fig.add_subplot(gs[0, 0])
    discrete_matshow(np.round(geo.G[:,:,0], 2), ax=ax0)

    # Temperature fields
    for idx, (temp, title) in enumerate([(T0, 'Reference Wall'), (T, geo.Nom_Geometry)]):
        ax = fig.add_subplot(gs[0, idx+1])
        
        # Convert temperatures to Celsius
        temp_C = np.flip(np.round(temp - 273.15, 1), axis=0)
        
        # Create temperature plot
        c = ax.contourf(temp_C, levels=np.linspace(temp_C.min(), temp_C.max(), 50), 
                       cmap='jet')
        fig.colorbar(c, ax=ax, label='Temperature (°C)')
        
        # Add boundary lines for wall and floor
        if idx == 1:  # Only for thermal bridge plot
            # Vertical wall line
            ax.axvline(x=p, color='white', linestyle='--', alpha=0.8, label='Wall boundary')
            
            # Horizontal floor lines
            ax.axhline(y=m-q1, color='red', linestyle='--', alpha=0.8, label='Upper floor boundary')
            ax.axhline(y=m-q2, color='red', linestyle='--', alpha=0.8, label='Lower floor boundary')
            
            # Add legend
            ax.legend(loc='upper right')
            
            # Add heat flux vectors
            for i in range(0, n, 10):
                for j in range(0, m, 10):
                    Jx, Jy = J[j, i]
                    scale = np.sqrt(Jx**2 + Jy**2)
                    if scale > 0:
                        ax.arrow(i, j, 5*Jx/scale, 5*Jy/scale, 
                                head_width=1, head_length=1, fc='white', ec='white',
                                alpha=0.5)
        
        # Set labels and title
        ax.set_title(f'Temperature field: {title}')
        ax.set_xlabel('x-axis (m)')
        ax.set_ylabel('y-axis (m)')
        
        # Set ticks
        ax.set_xticks(L_N[0::n//10]+L_N[-1:])
        ax.set_xticklabels([f'{x:.2f}' for x in np.arange(0, Lmax+dx, Lmax/10)])
        ax.set_yticks(L_M[0::n//10]+L_M[-1:])
        ax.set_yticklabels([f'{y:.2f}' for y in np.arange(0, Hmax+dy, Hmax/10)])

    # Print results
    print(f'Linear thermal bridge coefficient: {Psi:.4f} W/(m·K) (Geometry: {geo.Nom_Geometry})')
    print(f'Cost of thermal bridge losses: {prix*24*(Tint-Text2)*P*Psi:.4f} € per day')
    
    plt.show()

###################################### Commande à l'utilisateur #########################
print("Choose which activity you want to do:\n")
print("1: Simple simulation of a geometry\n")
print("2: Optimization of the insulation thickness on the floor\n")
print("3: Optimization of the position of the thermal bridge breaker or the insulation strip\n")
answer = input("Your choice (1, 2 or 3): ")

# Reference geometry: wall
geo0 = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, Hmax, P, 
                Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef)
geo0.mur()

################################################################### Application #####################################################################
if answer == '1':
    geo = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, Hmax, P, 
                Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef)
    geo.rupteur()  # Choix de la géométrie à simuler

    # T,Temps_stat = Initialisation(),0
    T0, Temps_stat = schema_jacobi(geo0,CL)
    T, Temps_stat = schema_jacobi(geo,CL)

    print('Time taken to reach steady state:', Temps_stat, 'min\n')

    J = Jth(geo,T)

    Psi, p, q1, q2 = Psi_lin(geo0,T0,geo,T,J)
    Psi = round(Psi,4)
    
    Affiche(geo,T0,T,J,Psi,p,q1,q2)

###################################### Réduction de Psi au maximum ###########################################
################################ Optimisation de l'épaisseur de l'isolant ################################
elif answer == '2':
    ''' Evolution of insulation thickness '''

    geo = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, Hmax, P, 
                Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef)
    method_name = 'planelle_isolant'

    li_min, li_max, n_steps = 0.05, 0.22, 5
    best_li, best_psi, L_li, L_psi = optimisation_floor(geo0, geo, li_min, li_max, n_steps, method_name)
    print(f"Optimal insulation thickness: {best_li:.4f} m")
    print(f"Corresponding linear thermal transmittance Ψ: {best_psi:.4f} W/m·K")

    plt.figure('Insulation Thickness Optimization')
    plt.plot(L_li, L_psi, marker='o')
    plt.title('Insulation Thickness Optimization')
    plt.xlabel('Insulation thickness (m)')
    plt.ylabel('Linear thermal transmittance Ψ (W/m·K)')
    plt.grid()
    plt.show()

elif answer == '3':
    ''' Position of insulation block or thermal break '''

    geo = Geometry(n, m, e, li, hi, h, dx, dy, dt, Text1, Text2, Tint, Lmax, Hmax, P, 
                Cth0, Cth1, Cth2, hm1, w, eps_iso, Coef)
    method_name = "rupteur"  # ou 'planelle_isolant'
    best_Coef, best_psi, L_Coef, L_psi = optimisation_planelle(geo0, geo, method_name)
    print(f"Optimal position of block/thermal break: {best_Coef:.2f} m")
    print(f"Corresponding linear thermal transmittance Ψ: {best_psi:.4f} W/(m·K)")

    plt.figure('Block/Thermal Break Position Optimization')
    plt.plot(L_Coef, L_psi, marker='o')
    plt.title('Block/Thermal Break Position Optimization')
    plt.xlabel('Block/Thermal break position (m)')
    plt.ylabel('Linear thermal transmittance Ψ (W/m·K)')
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
