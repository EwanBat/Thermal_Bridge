import numpy as np

class Geometrie:
    def __init__(self, n, m, e, li, h, hi, dx, dy, Cth0, Cth1, Cth2, Coef=0.5, ep_iso=0.05):
        self.n = n
        self.m = m
        self.e = e
        self.li = li
        self.h = h
        self.hi = hi
        self.dx = dx
        self.dy = dy
        self.Cth0 = Cth0
        self.Cth1 = Cth1
        self.Cth2 = Cth2
        self.Coef = Coef
        self.ep_iso = ep_iso
        self.Nom_Geometrie = ""
        self.G = None  # Attribut pour stocker la matrice

    def mur(self):
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
        self.Nom_Geometrie = 'murale inversée'
        self.G = np.zeros([self.m, self.n, 2])
        self.G[:, 0] = [[-1, self.Cth0] for _ in range(self.m)]
        for i in range(1, int(self.li/self.dx)+1):
            self.G[:, i] = [[2, self.Cth2] for _ in range(self.m)]
        for i in range(int(self.li/self.dx)+1, int((self.e+self.li)/self.dx)+1):
            self.G[:, i] = [[1, self.Cth1] for _ in range(self.m)]
        for i in range(int((self.e+self.li)/self.dx)+1, self.n):
            self.G[:, i] = [[0, self.Cth0] for _ in range(self.m)]

    def ss(self):
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