import unittest
import numpy as np

from PB import PB
from pbMatrices import pbMatrices
from Parabolic_System import Parabolic_System

theta = np.array([1.,.0046, 1.23])
tau = 5.
P = 11
k = 1
u = np.array(
    [0, 12.25, 24.5, 36.75, 49, 47, 45, 47.333, 49.667, 52, 49.25, 46.5, 43.75, 41, 39, 37, 36, 35, 34, 32.333, 30.667,
     29, 27.667, 26.333, 25, 24.333, 23.667, 23, 22.167, 21.333, 20.5, 19.667, 18.833, 18, 17.333, 16.667, 16, 15, 14,
     13, 12, 11, 10, 6.6667, 3.3333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
y = np.array(
    [0, 0, 0, 0, 1.7657, 6.1143, 11.289, 16.125, 19.03, 19.934, 20.982, 21.791, 24.744, 28.649, 32.077, 34.315,
     34.315, 32.839, 32.553, 32.22, 31.887, 31.41, 31.077, 30.744, 30.601, 30.41, 30.458, 30.077, 29.887, 29.649, 29.22,
     28.791, 28.553, 28.268, 28.125, 28.315, 28.172, 27.887, 27.649, 27.125, 26.934, 26.934, 26.649, 26.125, 25.41,
     24.22, 23.22, 22.125, 21.268, 20.553, 19.791, 19.315, 18.506, 17.601, 16.744, 16.077, 15.315, 14.934, 14.553,
     14.315, 13.839, 13.172, 12.363, 11.458, 10.934, 10.41, 10.03, 9.4105, 8.9343, 8.5533, 8.22, 7.8867, 7.5057, 6.9343,
     6.4581, 6.2676, 6.1248, 6.22, 6.3152, 6.1724, 5.8867, 5.5533, 5.2676, 5.1248, 5.3152, 5.1724, 4.8867, 4.5533,
     4.2676, 4.22, 4.1724, 3.8867, 3.5533, 3.2676, 3.22, 3.1724, 2.8867, 2.5533, 2.2676, 2.22, 2.1724, 1.8867, 1.5533,
     1.2676, 1.1248, 1.22, 1.22, 1.22, 1.3152, 1.1724, 0.88667, 0.55333, 0.26762, 0.12476, 0.22, 0.22, 0.22, 0.22, 0.22,
     0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.24095, 0.20952, 0.14667, 0.073333, 0.010476, 0, 0, 0, 0 ,0 ,0 ,0 ,0])

u = u[::10]
y = y[::10]

n = u.shape[0]
T = n*tau

Z = Parabolic_System(theta[1],theta[2], u, n, tau)
pbM = pbMatrices(theta, u, y, P, T, n, tau, k)
pb = PB()


class MyTest(unittest.TestCase):

    def test_plot_conv_kernel(self):
        import numpy as np
        import matplotlib.pyplot as plt
        x = np.arange(0., T, tau)
        y = np.array([pbM.conv_kernel(i) for i in x]).reshape(x.shape[0])
        plt.plot(x, y)
        plt.show()

    def test_vy_th_shape(self):
        self.assertEqual(pbM.vy_th.shape, (pbM.n, pbM.n))

    def test_vy_th_nonsingular(self):
        self.assertNotEqual(np.linalg.det(pbM.vy_th), 0.)

    def test_vy_th_nonsingular_diffthetas(self):
        pbM.set_theta_update_operators([1., .001, 0.1], u)
        self.assertNotEqual(np.linalg.det(pbM.vy_th), 0.)

    def test_vy_th_nonsingular_diffus(self):
        pbM.set_theta_update_operators(pbM.theta, np.ones(pb.n))
        self.assertNotEqual(np.linalg.det(pbM.vy_th), 0.)
        
    def test_vy_th_positive_semidefinite(self):
        pbM.set_theta_update_operators(pbM.theta, pbM.u)
        self.assertTrue(all([x>0. for x in np.linalg.eigvalsh(pbM.vy_th)]))

    def test_vy_th_real(self):
        from numpy import conj
        pbM.set_theta_update_operators(pbM.theta, pbM.u)
        self.assertTrue(all([x == conj(x) for x in pbM.vy_th.reshape(-1)]))

    def test_vy_th_symmetric(self):
        from numpy import conj
        #print np.linalg.eigvalsh(pbM.vy_th)
        self.assertTrue(all([x == conj(x) for x in np.linalg.eigvalsh(pbM.vy_th).reshape(-1)]))

    def test_vhat_positive_semidefinite(self):
        pbM.set_theta_update_operators(pbM.theta, pbM.u)
        print pbM.vhat
        print np.linalg.eigvalsh(pbM.vhat)
        self.assertTrue(all([x>0. for x in np.linalg.eigvalsh(pbM.vhat)]))

if __name__ == "__main__":
    unittest.main()
