import numpy as np

from geometrie import Geometry

def gradxy(geo: Geometry, T: np.ndarray, i: int, j: int) -> np.ndarray:
    """
    Calculates the temperature gradient at coordinates (i,j)
    
    Args:
        geo (Geometry): Object containing geometric parameters
        T (np.array): Temperature matrix
        i,j (int): Grid position indices
        
    Returns:
        np.array: Gradient [dT/dx, dT/dy]
    """
    if i == geo.n-1:  # Right edge
        if j == geo.m-1:  # Bottom right corner
            Y = (T[geo.m-2,i]-T[j,i])/geo.dy  # Backward finite difference in y
            X = (T[j,geo.n-2]-T[j,geo.n-1])/geo.dx  # Backward finite difference in x
        else:  # Right edge (not corner)
            Y = (T[j+1,i]-T[j,i])/geo.dy  # Forward finite difference in y
            X = (T[j,geo.n-2]-T[j,geo.n-1])/geo.dx  # Backward finite difference in x
    else:  # Interior points or left edge
        if j == geo.m-1:  # Bottom edge
            Y = (T[geo.m-2,i]-T[j,i])/geo.dy  # Backward finite difference in y
            X = (T[j,i+1]-T[j,i])/geo.dx  # Forward finite difference in x
        else:  # Interior points
            Y = (T[j+1,i]-T[j,i])/geo.dy  # Forward finite difference in y
            X = (T[j,i+1]-T[j,i])/geo.dx  # Forward finite difference in x
    return np.array([X,Y])

def Initialisation(geo: Geometry) -> np.ndarray:
    """
    Initializes temperature field with linear variation
    between Text1 (exterior) and Tint (interior)
    
    Args:
        geo (Geometry): Object containing geometric parameters
        
    Returns:
        np.array: Initial temperature matrix
    """
    T = np.zeros([geo.m, geo.n])
    for j in range(geo.m):
        T[j,:] = [geo.Text1 - (geo.Text1-geo.Tint)*geo.dx*i/geo.Lmax for i in range(geo.n)]
    return T

def calc(geo: Geometry, T: np.ndarray, i: int, j: int)-> float:
    """
    Calculates new temperature at point (i,j) using relaxation method
    
    Args:
        geo (Geometry): Object containing geometric parameters
        T (np.array): Current temperature matrix
        i,j (int): Grid position indices
        
    Returns:
        float: Calculated new temperature
    """
    # Relaxation method with coefficient w
    Term = geo.w*(1/4)*(T[j+1,i] + T[j-1,i] + T[j,i+1] + T[j,i-1]) + (1-geo.w)*T[j,i]
    return Term

def itere(geo: Geometry, T: np.ndarray, Condition) -> float:
    """
    Performs a complete iteration over the entire grid
    
    Args:
        geo (Geometry): Object containing geometric parameters
        T (np.array): Temperature matrix
        Condition (function): Function defining boundary conditions
        
    Returns:
        float: Sum of squared differences between two iterations
    """
    Sum = 0
    for i in range(geo.n):
        for j in range(geo.m):
            Bool, val = Condition(geo, T, i, j)  # Check if point is at boundary conditions
            if not Bool:  # Internal point: calculation by relaxation
                Term = calc(geo, T, i, j)
            else:  # Boundary point: use imposed value
                Term = val
            Sum += (Term - T[j,i])**2  # Accumulate difference
            T[j,i] = Term  # Update temperature
    return Sum

def schema_jacobi(geo: Geometry, Condition)-> tuple[np.ndarray, float]:
    """
    Solves heat equation using Jacobi method until convergence
    
    Args:
        geo (Geometry): Object containing geometric parameters
        Condition (function): Function defining boundary conditions
        
    Returns:
        tuple: (Final temperature matrix, Simulation time)
    """
    T = Initialisation(geo)  # Field initialization
    Ind = itere(geo,T,Condition)  # First iteration
    Counter = 1
    
    # Iterate until convergence (threshold 0.1)
    while Ind >= 1e-1:
        Ind = itere(geo,T,Condition)
        Counter += 1
        
    # Calculate simulation time and convert to Celsius
    Time_stat = round(Counter*geo.dt/60,1)
    T = T - 273
    
    return T,Time_stat