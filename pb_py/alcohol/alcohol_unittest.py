import unittest
import numpy as np


from PB import PB

pb = PB()
pb.calculate_new_operators(pb.thetas[:, 0], pb.us[:, 0])
pb.set_vhat(pb.thetas[:, 0])






class MyTest(unittest.TestCase):
    """
    def test_vy_th_shape(self):
        self.assertEqual(pb.vy_th.shape, (pb.n, pb.n))

    def test_lmatrix_nonzero(self):
        self.assertIsNot(pb.lmatrix, np.zeros(pb.lmatrix.shape))

    def test_vhat(self):
        print pb.vhat
    """

    def test_index(self):
        a = np.zeros((8, 50))
        a[3, 2] = 1
        self.assertEqual(a[3,2],a[[3,2]])

if __name__ == "__main__":
    unittest.main()
    from rkhs import Green1_eigen

    g = Green1_eigen(pb.P, pb.T, pb.thetas[:, 0])
    t = np.arange(0., pb.T, pb.tau)
    import matplotlib.pyplot as plt

    plt.plot(t, g.eifs[0](t))
    plt.plot(t, g.eifs[1](t))
    plt.plot(t, g.eifs[2](t))
    plt.show()