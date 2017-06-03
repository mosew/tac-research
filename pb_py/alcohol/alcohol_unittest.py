import unittest
import numpy as np

from PB import PB

pb = PB()
pb.calculate_new_operators(pb.thetas[:, 0], pb.us[:, 0])
pb.set_vhat(pb.thetas[:, 0])


class MyTest(unittest.TestCase):

    def test_vy_th_shape(self):
        self.assertEqual(pb.vy_th.shape, (pb.n, pb.n))

    def test_vy_th_nonsingular(self):
        self.assertNotEqual(np.linalg.det(pb.vy_th), 0.)

    def test_vy_th_nonsingular_diffthetas(self):
        pb.calculate_new_operators([10., .005, 0.5], pb.us[:, 0])
        self.assertNotEqual(np.linalg.det(pb.vy_th), 0.)

    def test_vy_th_nonsingular_diffus(self):
        pb.calculate_new_operators(pb.thetas[:, 0], np.ones(pb.n))
        self.assertNotEqual(np.linalg.det(pb.vy_th), 0.)

if __name__ == "__main__":
    unittest.main()
