import numpy as np

class Geometry:
    """
    Class representing different geometric configurations of thermal bridges.
    
    This class allows modeling several types of wall and insulation configurations
    for studying thermal bridges in buildings. It manages a 2D grid with
    different materials (air, concrete, insulation) characterized by their thermal properties.

    Attributes:
        n, m (int): Discretization grid dimensions
        e (float): Concrete wall thickness (m)
        li (float): Insulation thickness (m)
        h (float): Floor height (m)
        hi (float): Half-height of insulation (m)
        dx, dy (float): Spatial discretization steps
        dt (float): Time step
        Text1 (float): Initial external temperature (K)
        Text2 (float): Final external temperature (K)
        Tint (float): Interior temperature (K)
        Lmax (float): Maximum domain dimension (m)
        P (float): Room depth (m)
        Cth0 (float): Air thermal conductivity (W/m.K)
        Cth1 (float): Concrete thermal conductivity (W/m.K)
        Cth2 (float): Insulation thermal conductivity (W/m.K)
        hm1 (float): External surface heat transfer coefficient (W/m².K)
        w (float): Relaxation coefficient
        eps_iso (float): Floor insulation thickness (m)
        Coef (float): Block/thermal break position coefficient
        G (array): Thermal properties matrix [material_type, conductivity]
        Nom_Geometry (str): Name of geometric configuration
    """

    def __init__(self, n: int, m: int, 
                 e: float, li: float, h: float, hi: float, dx: float, dy: float, dt: float, 
                 Text1: float, Text2: float, Tint: float, Lmax: float, P: float,
                 Cth0: float, Cth1: float, Cth2: float, hm1: float,
                 w: float, eps_iso: float, Coef: float):
        """
        Initializes a new geometric configuration.
        
        Material types are coded as follows in matrix G:
            -1 : External air
             0 : Internal air
             1 : Concrete
             2 : Insulation
        """
        self.n = n
        self.m = m
        self.e = e
        self.li = li
        self.h = h
        self.hi = hi

        self.dx = dx
        self.dy = dy
        self.dt = dt

        self.Text1 = Text1
        self.Text2 = Text2
        self.Tint = Tint
        self.Lmax = Lmax
        self.P = P

        self.Cth0 = Cth0
        self.Cth1 = Cth1
        self.Cth2 = Cth2
        self.hm1 = hm1

        self.w = w
        self.eps_iso = eps_iso  # Épaisseur d'isolation dans le plancher (m)
        self.Coef = Coef
        
        self.Nom_Geometry = ""
        self.G = np.zeros((self.m,self.n,2), dtype=float)  # Attribut pour stocker la matrice

    ## Mise à jour des paramètres
    def update_params(self, **kwargs):
        """
        Updates geometry parameters.
        
        Allows modification of class attributes by passing key-value pairs.
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(f"{key} n'est pas un attribut valide de Geometry.")

    ## Configurations géométriques possibles

    def mur(self):
        """
        Creates a simple wall configuration with external insulation.
        Structure: [Ext | Concrete | Insulation | Int air]
        """
        self.Nom_Geometry = 'murale'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.e+self.li)/self.dx)+1):
            self.G[:, i] = [[2, self.Cth2] for _ in range(self.m)]
        for i in range(int((self.e+self.li)/self.dx)+1, self.n):
            self.G[:, i] = [[0, self.Cth0] for _ in range(self.m)]

    def mur_inv(self):
        """
        Creates a wall configuration with internal insulation.
        Structure: [Ext | Insulation | Concrete | Int air]
        """
        self.Nom_Geometry = 'murale inversée'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.li/self.dx)+1):
            self.G[:, i] = [[2, self.Cth2] for _ in range(self.m)]
        for i in range(int(self.li/self.dx)+1, int((self.e+self.li)/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int((self.e+self.li)/self.dx)+1, self.n):
            self.G[:, i] = [[0, self.Cth0] for _ in range(self.m)]

    def sans_isolation(self):
        """
        Creates a configuration without insulation.
        Structure: [Ext | Concrete | Int air]
        """
        self.Nom_Geometry = 'sans isolation'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int((self.e)/self.dx)+1, self.n):
            self.G[:int(self.hi/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]

    def classique(self):
        """
        Creates a classic configuration with distributed insulation.
        Structure: [Ext | Concrete | Insulation+Concrete | Int air]
        """
        self.Nom_Geometry = 'classique'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int(self.hi/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]

    def inv(self):
        """
        Creates an inverted configuration with internal insulation.
        Structure: [Ext | Insulation | Concrete+Floor | Int air]
        """
        self.Nom_Geometry = 'inversée'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.li/self.dx)+1):
            self.G[:, i] = [[2, self.Cth2] for _ in range(self.m)]
        for i in range(int(self.li/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int(self.hi/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.hi+self.h)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]

    def rupteur(self):
        """
        Creates a configuration with thermal break.
        Structure: [Ext | Concrete | Break | Concrete | Insulation+Concrete | Int air]
        The thermal break is an air zone between two concrete sections.
        """
        self.Nom_Geometry = 'rupteur thermique'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.Coef*self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.Coef*self.e/self.dx)+1, int((self.Coef+1/10)*self.e/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[1, self.Cth1] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.Coef+1/10)*self.e/self.dx)+1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int(self.hi/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]

    def planelle_isolant(self):
        """
        Creates a configuration with insulating block.
        Structure: [Ext | Concrete | Block | Concrete | Insulation+Concrete | Int air]
        The block is an insulation zone between two concrete sections.
        """
        self.Nom_Geometry = 'planelle isolant'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.Coef*self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.Coef*self.e/self.dx)+1, int((self.Coef+1/10)*self.e/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[1, self.Cth1] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.Coef+1/10)*self.e/self.dx)+1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int(self.hi/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]

    def classique_plancher_isolant(self):
        """
        Creates a configuration with insulated floor.
        Structure: [Ext | Concrete | Insulation+Concrete | Air+Insulation+Concrete+Insulation | Int air]
        Includes an insulation layer in the floor.
        """
        self.Nom_Geometry = 'plancher isolé'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int((self.hi-self.eps_iso)/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int((self.hi-self.eps_iso)/self.dy)+1)]
            self.G[int((self.hi-self.eps_iso)/self.dy)+1:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int((self.hi-self.eps_iso)/self.dy)+1, int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:int((self.h+self.hi+self.eps_iso)/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, int((self.h+self.hi+self.eps_iso)/self.dy)+1)]
            self.G[int((self.h+self.hi+self.eps_iso)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi+self.eps_iso)/self.dy)+1, self.m)]

    def classique_plancher_simple(self):
        """
        Creates a configuration with simple floor.
        Structure: [Ext | Concrete | Insulation+Concrete | Air+Insulation+Concrete | Int air]
        Simplified version without additional floor insulation.
        """
        self.Nom_Geometry = 'plancher simple'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int((self.hi-self.eps_iso)/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int((self.hi-self.eps_iso)/self.dy)+1)]
            self.G[int((self.hi-self.eps_iso)/self.dy)+1:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int((self.hi-self.eps_iso)/self.dy)+1, int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]

