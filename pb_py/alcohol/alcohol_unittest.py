import unittest
import numpy as np

from PB import PB
from pbMatrices import pbMatrices
from Parabolic_System import Parabolic_System

from scipy.stats import exponnorm


theta = np.array([1., .0046, 1.23])
tau = 1./12.
n = 46
T = n*tau
P = 7
k = 55
u = np.array(
    [.01 * exponnorm.pdf((x - 30.) / 40., 0.3) for x in np.arange(tau, (n + 1) * tau, tau)])
y = np.array(
    [.01 * exponnorm.pdf((x - 30.) / 40., 0.5) for x in np.arange(tau, (n + 1) * tau, tau)])

u_final = np.array(
    [400. * exponnorm.pdf((x - 30.) / 4https://giphy.com/gifs/C0mVQRWkvO7kI/html5., 0.3) for x in np.linspace(tau, n * tau, n)])
y_final = np.ones(u_final.shape)

Z = Parabolic_System(theta[0],theta[1],u, n)
pbM = pbMatrices(theta, u, y, P, T, n, tau, k)
pb = PB()



class MyTest(unittest.TestCase):

    def test_plot(self):
        from scipy.stats import exponnorm
        import matplotlib.pyplot as plt
        K=1.5
        x = np.linspace(exponnorm.ppf(0.01, K), exponnorm.ppf(0.99, K), 100)
        plt.plot(x, exponnorm.pdf(x, K))

    def test_vy_th_shape(self):
        self.assertEqual(pbM.vy_th.shape, (pbM.n, pbM.n))

    def test_vy_th_nonsingular(self):
        self.assertNotEqual(np.linalg.det(pbM.vy_th), 0.)

    def test_vy_th_nonsingular_diffthetas(self):
        pbM.set_theta_update_operators([1., .001, 0.1],u,y)
        self.assertNotEqual(np.linalg.det(pbM.vy_th), 0.)

    def test_vy_th_nonsingular_diffus(self):
        pbM.set_theta_update_operators(pbM.theta, np.ones(pb.n), y)
        self.assertNotEqual(np.linalg.det(pbM.vy_th), 0.)
        
    def test_vy_th_positive_semidefinite(self):
        pbM.set_theta_update_operators(pbM.theta, pbM.u, pbM.y)
        self.assertTrue(all([x>0. for x in np.linalg.eigvalsh(pbM.vy_th)]))

    def test_vy_th_real(self):
        from numpy import conj
        pbM.set_theta_update_operators(pbM.theta, pbM.u, pbM.y)
        self.assertTrue(all([x == conj(x) for x in pbM.vy_th.reshape(-1)]))

    def test_vy_th_symmetric(self):
        from numpy import conj
        print np.linalg.eigvalsh(pbM.vy_th)
        self.assertTrue(all([x == conj(x) for x in np.linalg.eigvalsh(pbM.vy_th).reshape(-1)]))

    def test_vhat_positive_semidefinite(self):
        pbM.set_theta_update_operators(pbM.theta, pbM.u, pbM.y)
        self.assertTrue(all([x>0. for x in np.linalg.eigvalsh(pbM.vhat)]))

if __name__ == "__main__":
    unittest.main()
