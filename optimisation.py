import numpy as np

from geometrie import Geometry
from fonction_calcul import schema_jacobi  # si dans le même package
from fonction_therm import CL, Jth, Psi_lin

def optimisation_floor(geo0: Geometry, geo: Geometry, e_min: float, e_max: float, n_steps: int) -> tuple:
    """
    Optimizes the thickness of the insulation in a floor thermal bridge configuration to minimize energy losses.

    Parameters:
    - geom: Geometrie object containing the geometry and material properties.
    - e_min: Minimum thickness of the insulation (in meters).
    - e_max: Maximum thickness of the insulation (in meters).
    - n_steps: Number of steps to divide the thickness range into.

    Returns:
    - A tuple containing:
        - Optimal thickness of the insulation (in meters).
        - Corresponding thermal bridge coefficient Ψ (in W/m·K).
        - Estimated annual energy loss (in kWh/m).
    """
    best_e = float(e_min)
    best_psi = float('inf')
    L_e = np.linspace(e_min, e_max, n_steps)
    L_psi = []

    T0 = schema_jacobi(geo0, CL)

    for e in L_e:
        geo.update_params(e=e)
        T = schema_jacobi(geo, CL)
        J_th = Jth(T, geo)
        psi = Psi_lin(geo0, T0, geo, T, J_th)
        L_psi.append(psi)
        

        if psi < best_psi:
            best_psi = psi
            best_e = e

    return best_e, best_psi, np.array(L_e), np.array(L_psi)