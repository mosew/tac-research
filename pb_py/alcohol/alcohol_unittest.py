import unittest
import numpy as np
from PB import *

pb = PB()


class MyTest(unittest.TestCase):
    def test_compute_li(self):
        print pb.compute_L_i([1., .0046, 1.23], pb.u, 3)

if __name__ == "__main__":
    unittest.main()
