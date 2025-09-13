import numpy as np

from geometrie import Geometry
from fonction_calcul import schema_jacobi  # si dans le même package
from fonction_therm import CL, Jth, Psi_lin

def optimisation_floor(geo0: Geometry, geo: Geometry, li_min: float, li_max: float, n_steps: int, method_name: str) -> tuple:
    """
    Optimizes the thickness of the insulation in a floor thermal bridge configuration to minimize energy losses.

    Parameters:
    - geom: Geometrie object containing the geometry and material properties.
    - li_min: Minimum thickness of the insulation (in meters).
    - li_max: Maximum thickness of the insulation (in meters).
    - n_steps: Number of steps to divide the thickness range into.

    Returns:
    - A tuple containing:
        - Optimal thickness of the insulation (in meters).
        - Corresponding thermal bridge coefficient Ψ (in W/m·K).
        - Array of thickness values tested (in meters).
        - Array of corresponding Ψ values (in W/m·K).
    """
    best_li = float(li_min)
    best_psi = float('inf')
    L_li = np.linspace(li_min, li_max, n_steps)
    L_psi = []

    T0, Time_stat = schema_jacobi(geo0, CL)

    for li in L_li:
        geo.update_params(li=li)
        geo.update_geometry(method_name) # Update geometry with new insulation thickness

        T, Time_stat = schema_jacobi(geo, CL)
        J_th = Jth(geo, T)
        psi, p, q1, q2 = Psi_lin(geo0, T0, geo, T, J_th)
        L_psi.append(psi)
        
        if psi < best_psi:
            best_psi = psi
            best_li = li
        
        print("Optimisation at ", round(100 * (li - li_min) / (li_max - li_min),1), "% completed", end="\r")

    return best_li, best_psi, L_li, np.array(L_psi)

def optimisation_planelle(geo0: Geometry, geo: Geometry, method_name: str)-> tuple:
    """
    Optimizes the position of the insulation in a wall thermal bridge configuration to minimize energy losses.

    Parameters:
    - geom: Geometrie object containing the geometry and material properties.

    Returns:
    - A tuple containing:
        - Optimal position coefficient of the insulation (between 0 and 1).
        - Corresponding thermal bridge coefficient Ψ (in W/m·K).
        - Array of position coefficients tested (between 0 and 1).
        - Array of corresponding Ψ values (in W/m·K).
    """
    best_Coef = 0.01
    best_psi = float('inf')
    L_Coef = np.arange(0.1, 1, 0.2)
    L_psi = []
    T0, Time_stat = schema_jacobi(geo0, CL)

    for Coef in L_Coef:
        geo.update_params(Coef=Coef)
        geo.update_geometry(method_name) # Update geometry with new position coefficient

        T, Time_stat = schema_jacobi(geo, CL)
        J_th = Jth(geo, T)
        psi, p, q1, q2 = Psi_lin(geo0, T0, geo, T, J_th)
        L_psi.append(psi)
        
        if psi < best_psi:
            best_psi = psi
            best_Coef = Coef
        
        print("Optimisation at ", round(100 * (Coef - 0.1) / (0.9 - 0.1),1), "% completed", end="\r")
    
    print(L_psi)
    return best_Coef, best_psi, L_Coef, np.array(L_psi)