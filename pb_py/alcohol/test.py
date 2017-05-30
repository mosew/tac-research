import unittest

from scipy.stats import beta, norm
import numpy as np

class PBTestCase(unittest.TestCase):
    def test_compute_L_i(self):
        f = lambda x: norm.rvs(0.3 * beta.pdf(x, 12, 7) + 0.6 * beta.pdf(x, 4, 11), scale=0.05)
        i = 3

    def slice(unittest.TestCase):
