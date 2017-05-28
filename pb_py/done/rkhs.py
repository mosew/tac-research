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
        self.eivs = list(P)
        self.eifs = list(P)

        from math import pi, sqrt, sin
        from sympy import symbols, lambdify

        s, theta = symbols("s,theta")
        for j in range(1, self.P + 1):
            self.eivs[j] = self.theta[0] * (self.T ** 2) / (((j - 1) * pi) + pi / 2) ** 2
            self.eifs[j] = lambdify(s, sqrt(2 / self.T) * sin((s / self.T) * (j * pi - pi / 2)), "numpy")

    def set_theta(self,theta):
        for j in range(1,self.P + 1):
            self.eivs[j] = theta[0] * (self.T ** 2) / (((j - 1) * pi) + pi / 2) ** 2