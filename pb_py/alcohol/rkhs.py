# rkhs


class Green1_eigen(object):
    """
    Returns an object with two useful attributes:
        eiv: a 1xP numpy array containing the P largest eigenvalues of the reproducing kernel
        eif: a 1xP list containing the eigenfunctions corresponding to these largest eigenvalues

    """

    def __init__(self, P, T, theta):
        self.P = P
        self.T = T
        self.theta = theta
        self.eifs = [0.0] * P

        import numpy as np
        self.eivs = np.zeros(P)
        self.deivs = np.zeros((P, theta.shape[0]))
        self.d2eivs = np.zeros((P, theta.shape[0], theta.shape[0]))

        from math import pi, sqrt, sin

        for j in range(self.P):
            self.eivs[j] = self.theta[0] * (self.T ** 2) / ((j * pi) + pi / 2) ** 2
            self.eifs[j] = lambda s: sqrt(2 / self.T) * sin((s / self.T) * (j * pi - pi / 2))
            self.deivs[j] = (self.T ** 2) / ((j * pi) + pi / 2) ** 2
            self.d2eivs[j] = 0

    def set_theta(self,theta):
        from math import pi
        for j in range(self.P):
            self.eivs[j] = theta[0] * (self.T ** 2) / ((j * pi) + pi / 2) ** 2
            self.deivs[j] = (self.T ** 2) / ((j * pi) + pi / 2) ** 2
            self.d2eivs[j] = 0