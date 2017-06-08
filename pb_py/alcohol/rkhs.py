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
        self.eifs = [lambda x: 0] * P

        import numpy as np
        self.eivs = np.zeros(P)
        self.deivs = np.zeros((P, theta.shape[0]))
        self.d2eivs = np.zeros((P, theta.shape[0], theta.shape[0]))

        from numpy import pi, sqrt, sin

        self.eifs = [(lambda s, j=j: (sqrt(2./self.T) * sin((s/self.T) * ((j+1) * pi - pi / 2.)))) for j in range(self.P)]

        for j in np.arange(1, self.P + 1):
            self.eivs[j-1] = 1.
            # self.eivs[j-1] = self.theta[0] * (self.T ** 2) / ((j * pi) + pi / 2) ** 2
            # self.deivs[j-1, 0] = (self.T ** 2) / ((j * pi) + pi / 2) ** 2

    def set_theta(self,theta):
        from numpy import pi
        for j in range(self.P):
            self.eivs[j-1]= 1.
            # self.eivs[j] = theta[0] * (self.T ** 2) / ((j * pi) + pi / 2) ** 2
            # self.deivs[j] = (self.T ** 2) / ((j * pi) + pi / 2) ** 2
            # self.d2eivs[j] = 0.


if __name__ == "__main__":
    import numpy as np
    g = Green1_eigen(8, 100., np.array([1., 1., 1.]))
    t = np.arange(0., 100., 5.)
    import matplotlib.pyplot as plt

    plt.plot(t, g.eifs[0](t))
    plt.plot(t, g.eifs[1](t))
    plt.plot(t, g.eifs[2](t))
    plt.show()
