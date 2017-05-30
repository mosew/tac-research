from PB import *
from scipy.stats import beta, norm
import numpy as np

pb = PB()


# Test compute_L_i
theta = [1., 10.]
f = lambda x:  norm.rvs(0.3 * beta.pdf(x, 12, 7) + 0.6 * beta.pdf(x, 4, 11), scale = 0.05)
i = 3

pb.calculate_new_operators(pb.thetas[:, 0])
pb.set_vhat(pb.thetas[:, 0])

c = np.linalg.cholesky(pb.vhat)
print c.shape
print c
