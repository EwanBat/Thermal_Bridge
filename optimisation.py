import numpy as np

from geometrie import Geometry
from fonction_calcul import schema_jacobi  # si dans le même package
from fonction_therm import CL, Jth, Psi_lin

def optimisation_floor(geo0: Geometry, geo: Geometry, li_min: float, li_max: float, n_steps: int) -> tuple:
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
        - Estimated annual energy loss (in kWh/m).
    """
    best_li = float(li_min)
    best_psi = float('inf')
    L_li = np.linspace(li_min, li_max, n_steps)
    L_psi = []

    T0, Time_stat = schema_jacobi(geo0, CL)

    for li in L_li:
        geo.update_params(li=li)
        T, Time_stat = schema_jacobi(geo, CL)
        J_th = Jth(geo, T)
        psi, p, q1, q2 = Psi_lin(geo0, T0, geo, T, J_th)
        L_psi.append(psi)
        
        if psi < best_psi:
            best_psi = psi
            best_li = li

    return best_li, best_psi, np.array(L_li), np.array(L_psi)