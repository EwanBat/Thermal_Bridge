import numpy as np

class Geometrie:
    """
    Classe représentant différentes configurations géométriques de ponts thermiques.
    
    Cette classe permet de modéliser plusieurs types de configurations de murs et d'isolation
    pour l'étude des ponts thermiques dans le bâtiment. Elle gère une grille 2D avec 
    différents matériaux (air, béton, isolant) caractérisés par leurs propriétés thermiques.

    Attributes:
        n, m (int): Dimensions de la grille de discrétisation
        e (float): Épaisseur du mur en béton (m)
        li (float): Épaisseur de l'isolant (m)
        h (float): Hauteur du plancher (m)
        hi (float): Demi-hauteur de l'isolant (m)
        dx, dy (float): Pas de discrétisation spatial
        dt (float): Pas de temps
        Text1 (float): Température extérieure initiale (K)
        Text2 (float): Température extérieure finale (K) 
        Tint (float): Température intérieure (K)
        Lmax (float): Dimension maximale du domaine (m)
        P (float): Profondeur de la pièce (m)
        Cth0 (float): Conductivité thermique de l'air (W/m.K)
        Cth1 (float): Conductivité thermique du béton (W/m.K)
        Cth2 (float): Conductivité thermique de l'isolant (W/m.K)
        hm1 (float): Coefficient d'échange surfacique extérieur (W/m².K)
        w (float): Coefficient de relaxation
        eps_iso (float): Épaisseur d'isolation dans le plancher (m)
        Coef (float): Coefficient de position planelle/rupteur
        G (array): Matrice des propriétés thermiques [type_matériau, conductivité]
        Nom_Geometrie (str): Nom de la configuration géométrique
    """

    def __init__(self, n, m, e, li, h, hi, dx, dy, dt, 
                 Text1, Text2, Tint, Lmax, P,
                 Cth0, Cth1, Cth2, hm1,
                 w, eps_iso, Coef):
        """
        Initialise une nouvelle configuration géométrique.
        
        Les types de matériaux sont codés comme suit dans la matrice G:
            -1 : Air extérieur
             0 : Air intérieur 
             1 : Béton
             2 : Isolant
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
        
        self.Nom_Geometrie = ""
        self.G = None  # Attribut pour stocker la matrice

    ## Mise à jour des paramètres
    def update_params(self, **kwargs):
        """
        Met à jour les paramètres de la géométrie.
        
        Permet de modifier les attributs de la classe en passant des paires clé-valeur.
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(f"{key} n'est pas un attribut valide de Geometrie.")

    ## Configurations géométriques possibles

    def mur(self):
        """
        Crée une configuration de mur simple avec isolation extérieure.
        Structure: [Ext | Béton | Isolant | Air int]
        """
        self.Nom_Geometrie = 'murale'
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
        Crée une configuration de mur avec isolation côté intérieur.
        Structure: [Ext | Isolant | Béton | Air int]
        """
        self.Nom_Geometrie = 'murale inversée'
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
        Crée une configuration sans isolation.
        Structure: [Ext | Béton | Air int]
        """
        self.Nom_Geometrie = 'sans isolation'
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
        Crée une configuration classique avec isolation répartie.
        Structure: [Ext | Béton | Isolant+Béton | Air int]
        """
        self.Nom_Geometrie = 'classique'
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
        Crée une configuration inversée avec isolation côté intérieur.
        Structure: [Ext | Isolant | Béton+Plancher | Air int]
        """
        self.Nom_Geometrie = 'inversée'
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
        Crée une configuration avec rupteur thermique.
        Structure: [Ext | Béton | Rupteur | Béton | Isolant+Béton | Air int]
        Le rupteur est une zone d'air entre deux sections de béton.
        """
        self.Nom_Geometrie = 'rupteur thermique'
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
        Crée une configuration avec planelle isolante.
        Structure: [Ext | Béton | Planelle | Béton | Isolant+Béton | Air int]
        La planelle est une zone d'isolant entre deux sections de béton.
        """
        self.Nom_Geometrie = 'planelle isolant'
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
        Crée une configuration avec plancher isolé.
        Structure: [Ext | Béton | Isolant+Béton | Air+Isolant+Béton+Isolant | Air int]
        Inclut une couche d'isolation dans le plancher.
        """
        self.Nom_Geometrie = 'plancher isolé'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int((self.hi-self.ep_iso)/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int((self.hi-self.ep_iso)/self.dy)+1)]
            self.G[int((self.hi-self.ep_iso)/self.dy)+1:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int((self.hi-self.ep_iso)/self.dy)+1, int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:int((self.h+self.hi+self.ep_iso)/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, int((self.h+self.hi+self.ep_iso)/self.dy)+1)]
            self.G[int((self.h+self.hi+self.ep_iso)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi+self.ep_iso)/self.dy)+1, self.m)]

    def classique_plancher_simple(self):
        """
        Crée une configuration avec plancher simple.
        Structure: [Ext | Béton | Isolant+Béton | Air+Isolant+Béton | Air int]
        Version simplifiée sans isolation supplémentaire dans le plancher.
        """
        self.Nom_Geometrie = 'plancher simple'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.e/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int(self.e/self.dx)+1, int((self.li+self.e)/self.dx)+1):
            self.G[:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[2, self.Cth2] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]
        for i in range(int((self.li+self.e)/self.dx)+1, self.n):
            self.G[:int((self.hi-self.ep_iso)/self.dy)+1, i] = [[0, self.Cth0] for _ in range(int((self.hi-self.ep_iso)/self.dy)+1)]
            self.G[int((self.hi-self.ep_iso)/self.dy)+1:int(self.hi/self.dy)+1, i] = [[2, self.Cth2] for _ in range(int((self.hi-self.ep_iso)/self.dy)+1, int(self.hi/self.dy)+1)]
            self.G[int(self.hi/self.dy)+1:int((self.h+self.hi)/self.dy)+1, i] = [[1, self.Cth1] for _ in range(int(self.hi/self.dy)+1, int((self.h+self.hi)/self.dy)+1)]
            self.G[int((self.h+self.hi)/self.dy)+1:, i] = [[0, self.Cth0] for _ in range(int((self.h+self.hi)/self.dy)+1, self.m)]

    