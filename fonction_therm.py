import numpy as np
from fonction_calcul import gradxy
from geometrie import Geometry


################################################ Thermodynamic Functions ###################################################
def CL(geo, T, i, j):
    """
    Determines boundary conditions for a point (j, i) in the grid
    based on geometry 'geo' and physical conditions.

    Parameters
    ----------
    geo : Geometry
        Contains geometric grid and material properties.
        Access to materials via geo.G[y, x] -> [type, conductivity]
    T : 2D array
        Temperature array.
    i : int
        x coordinate (column index).
    j : int
        y coordinate (row index).

    Returns
    -------
    tuple (bool, float)
        - bool : True if a boundary condition is applied, False otherwise.
        - float : Calculated temperature or 0 if no condition.
    """

    # --- LEFT EDGE: Fixed external temperature ---
    if i == 0:
        Tij = geo.Text2
        return True, Tij

    # --- RIGHT EDGE: Symmetric conditions or interior temperature ---
    elif i == geo.n - 1:
        if j == 0:  # Upper right corner
            if geo.G[j, i][0] != 0:  # Non-air → symmetry
                Tij = geo.w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j + 1, i] + T[geo.m - 2, i]) \
                      + (1 - geo.w) * T[j, i]
            else:  # Air → interior temperature
                Tij = geo.Tint

        elif j == geo.m - 1:  # Lower right corner
            if geo.G[j, i][0] != 0:  # Non-air → symmetry
                Tij = geo.w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j - 1, i] + T[1, i]) \
                      + (1 - geo.w) * T[j, i]
            else:  # Air → interior temperature
                Tij = geo.Tint

        else:  # Right edge without being in a corner
            if geo.G[j, i][0] != 0:
                Tij = geo.w * (1 / 4) * (T[j, i - 1] + T[j, i - 1] + T[j + 1, i] + T[j - 1, i]) \
                      + (1 - geo.w) * T[j, i]
            else:
                Tij = geo.Tint

        return True, Tij

    # --- UPPER EDGE ---
    elif j == 0:
        Tij = geo.w * (1 / 4) * (T[j, i + 1] + T[j, i - 1] + T[j + 1, i] + T[geo.m - 2, i]) \
              + (1 - geo.w) * T[j, i]
        return True, Tij

    # --- LOWER EDGE ---
    elif j == geo.m - 1:
        Tij = geo.w * (1 / 4) * (T[j, i + 1] + T[j, i - 1] + T[j - 1, i] + T[1, i]) \
              + (1 - geo.w) * T[j, i]
        return True, Tij

    # --- INTERFACE BETWEEN TWO MATERIALS ON X AXIS ---
    elif geo.G[j, i][0] != geo.G[j, i + 1][0] or geo.G[j, i][0] != geo.G[j, i - 1][0]:
        lamb1 = geo.G[j, i - 1][1]
        lamb2 = geo.G[j, i + 1][1]

        # Case with convection (type -1)
        if geo.G[j, i - 1][0] == -1 and geo.G[j, i][0] != geo.G[j, i - 1][0]:
            Tij = (lamb2 * T[j, i + 1] / geo.dx + T[j, i - 1] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        elif geo.G[j, i + 1][0] == -1 and geo.G[j, i][0] != geo.G[j, i + 1][0]:
            Tij = (lamb2 * T[j, i - 1] / geo.dx + T[j, i + 1] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        else:  # Purely conductive interface
            Tij = (lamb1 * T[j, i - 1] + lamb2 * T[j, i + 1]) / (lamb2 + lamb1)

        return True, Tij

    # --- INTERFACE BETWEEN TWO MATERIALS ON Y AXIS ---
    elif geo.G[j, i][0] != geo.G[j + 1, i][0] or geo.G[j, i][0] != geo.G[j - 1, i][0]:
        lamb1 = geo.G[j - 1, i][1]
        lamb2 = geo.G[j + 1, i][1]

        # Case with convection (type -1)
        if geo.G[j - 1, i][0] == -1 and geo.G[j, i][0] != geo.G[j - 1, i][0]:
            Tij = (lamb2 * T[j + 1, i] / geo.dx + T[j - 1, i] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        elif geo.G[j + 1, i][0] == -1 and geo.G[j, i][0] != geo.G[j + 1, i][0]:
            Tij = (lamb2 * T[j - 1, i] / geo.dx + T[j + 1, i] * geo.hm1) / (lamb2 / geo.dx + geo.hm1)
        else:  # Purely conductive interface
            Tij = (lamb1 * T[j - 1, i] + lamb2 * T[j + 1, i]) / (lamb2 + lamb1)

        return True, Tij

    # --- DEFAULT CASE: no special treatment, continue iteration ---
    else:
        return False, 0

def Jth(geo, T):
    """
    Calculates the thermal flux field J in the grid,
    according to Fourier's law:
        J = -λ * ∇T

    Parameters
    ----------
    geo : Geometry
        Geometry containing material properties.
        Access via geo.G[y, x, :] -> [type, conductivity]
    T : 2D array
        Temperature array.

    Returns
    -------
    J : 3D numpy array
        Thermal flux field of size (m, n, 2)
        - J[:, :, 0] = y component
        - J[:, :, 1] = x component
    """
    # Initialize flux array:
    # dimensions = [height (m), width (n), direction (2)]
    J = np.zeros((geo.m, geo.n, 2))

    # Loop through each grid point
    for i in range(geo.n):       # x axis
        for j in range(geo.m):   # y axis
            # Local temperature gradient (∇T)
            Grad = gradxy(geo, T, i, j)

            # Apply Fourier's law
            J[j, i] = -geo.G[j, i, 1] * Grad

    return J

def Psi_lin(geo0: Geometry, T0: np.ndarray, geo: Geometry, T: np.ndarray, J: np.ndarray) -> tuple[float, int, int, int]:
    """
    Calculates the linear thermal transmittance coefficient Psi to compare heat flows
    between a "reference" wall and the thermal bridge under study.

    Parameters
    ----------
    geo0 : Geometry
        Reference geometry (wall without thermal bridge).
        Access to materials via geo0.G[y, x, :] -> [type, conductivity]
    T0 : 2D array
        Temperature array for wall-only case.
    geo : Geometry
        Geometry of the model with thermal bridge.
    T : 2D array
        Temperature array for thermal bridge case.
    J : 3D array
        Calculated heat flux array (J[y, x, direction]).
        direction = 0 for y, 1 for x.

    Returns
    -------
    float
        Linear thermal transmittance coefficient Psi (W/m·K).
    """
    # --- Initialization ---
    i = 1  # Start after exterior air
    L_U0j = []        # Conductive flux per wall segment
    L_deltaT = []     # Temperature differences per segment
    L_Longueur = []   # Associated lengths

    compteur = 1      # Cumulative thickness of traversed material
    T_0 = T0[0, i]    # Temperature on exterior side
    T_lim = T0[0, i]  # Initial limit for temperature differences

    # --- Calculate reference wall linear fluxes ---
    while i < geo.n - 1:  # Until heated interior air
        # If material changes, calculate flux for previous segment
        if geo0.G[1, i + 1, 0] != geo0.G[1, i, 0]:
            L0i = compteur * geo.dx

            # Add linear flux (electrical analogy)
            L_U0j.append(geo0.G[1, i, 1] / L0i)
            L_deltaT.append(T0[0, i] - T_lim)
            L_Longueur.append(L0i)

            # Reset for new material
            compteur = 1
            T_lim = T0[0, i]
        else:
            compteur += 1
        i += 1

    # Convert to numpy arrays for vectorized calculations
    L_U0j = np.array(L_U0j)
    L_deltaT = np.array(L_deltaT)
    L_Longueur = np.array(L_Longueur)

    # --- Detection of characteristic indices ---
    if geo.Nom_Geometry == 'plancher simple':
        p = int((geo.li + geo.e) / geo.dx + 1)
        q1 = int((geo.hi - geo.eps_iso) / geo.dy)
        q2 = int((geo.hi + geo.h) / geo.dy + 1)
    else:
        p = int((geo.li + geo.e) / geo.dx + 1)
        q1 = int((geo.hi - geo.eps_iso) / geo.dy)
        q2 = int((geo.hi + geo.h + geo.eps_iso) / geo.dy + 1)

    T_lim = T0[0, i]  # Interior side temperature

    # --- Calculate linear fluxes ---
    Phi = geo.P * (
        np.sum(J[:q1, p, 0]) * geo.dy / geo.hi +
        np.sum(J[q2:, p, 0]) * geo.dy / geo.hi -
        np.sum(np.abs(J[q1 + 1, p:, 1])) * geo.dx / geo.Lmax -
        np.sum(np.abs(J[q2 - 1, p:, 1])) * geo.dx / geo.Lmax
    )

    # Reference flux (wall only)
    Phi_default = -np.sum(L_U0j * L_deltaT) * geo.P

    # --- Calculate Psi ---
    Psi = -(Phi - Phi_default) / (T_lim - T_0)

    return Psi, p, q1, q2

