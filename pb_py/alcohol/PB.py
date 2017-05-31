class PB(object):

    def __init__(self):

        self.P = 8
        self.T = 5*276
        self.n = 276
        self.tau = 5.
        self.burnin = 4
        self.K = 10
        self.alph = 3333
        self.nTheta = 3
        # What step of the MCMC are we on?
        self.k = 1

        import numpy as np
        self.infoTheta = np.array([[0, 0, 0],[1.6e-7, 0, 0],[0, 1e-2, 0]])
        self.thetas = np.zeros((self.nTheta, self.K))
        self.thetas[:, 0] = [1., 1., 1.]
        self.t = np.linspace(self.tau, self.T, self.n)
        self.a = np.zeros((self.P, self.K - self.burnin))
        self.ys = np.zeros((self.n, self.K - self.burnin))
        self.y = np.zeros(self.n)

        # Important matrices to be calculated
        self.vhat = np.zeros((self.nTheta, self.nTheta))
        self.lmatrix = np.zeros((self.P, self.n))
        self.dlmatrix = np.zeros((self.P, self.n, self.nTheta))
        self.d2lmatrix = np.zeros((self.P, self.n, self.nTheta, self.nTheta))
        self.vy_th = np.zeros((self.n, self.n))
        self.dvy_th = np.zeros((self.n, self.n, self.nTheta))
        self.d2vy_th = np.zeros((self.n, self.n, self.nTheta, self.nTheta))
        self.vy_thu = np.zeros((self.n, self.n))
        self.dvy_thu = np.zeros((self.n, self.n, self.nTheta))
        self.d2vy_thu = np.zeros((self.n, self.n, self.nTheta, self.nTheta))
        self.d2logpy_th = np.zeros((self.nTheta, self.nTheta))

        from rkhs import Green1_eigen
        self.g = Green1_eigen(self.P,self.T,self.thetas[:, 0])

        # This is used for the initial calculation of V(y|theta, u) and its derivatives wrt theta. It should be such that the operators are nonsingular.
        from scipy.stats import exponnorm
        self.u = lambda x: 100.*exponnorm.pdf(x, 0.5)

        from Parabolic_System import Parabolic_System
        self.Z = Parabolic_System(self.thetas[1,0],self.thetas[2,0])
        self.Z.define_operators()

        self.conv_kernel = lambda x: np.dot(self.Z.CNhat, np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau))), self.Z.BNhat))
        self.dconv_kernel_th1 = lambda x: np.dot(self.Z.CNhat, int(np.floor_divide(x, self.tau)) * np.dot(self.Z.dANhat_dq1, np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau)-1)), self.Z.BNhat)))
        self.d2conv_kernel_th1 = lambda x: np.dot(self.Z.CNhat,
                                                  np.dot(int(np.floor_divide(x, self.tau)) *
                                                        (np.dot(self.Z.d2ANhat_d2q1,
                                                                np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x,self.tau))-1)) +
                                                                (np.floor_divide(x, self.tau) - 1) * np.dot(self.Z.dANhat_dq1,
                                                                                                            np.linalg.matrix_power(self.Z.ANhat,
                                                                                                                                   int(np.floor_divide(x, self.tau)) - 2)
                                                                                                            )
                                                        ),
                                                        self.Z.BNhat)
                                                 )

    def recalculate_kernel(self,q1,q2):
        from Parabolic_System import Parabolic_System
        self.Z = Parabolic_System(q1,q2)
        self.Z.define_operators()

        self.conv_kernel = lambda x: np.dot(self.Z.CNhat, np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau))), self.Z.BNhat))
        self.dconv_kernel_th1 = lambda x: np.dot(self.Z.CNhat, int(np.floor_divide(x, self.tau)) * np.dot(self.Z.dANhat_dq1, np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau)-1)), self.Z.BNhat)))
        self.d2conv_kernel_th1 = lambda x: np.dot(self.Z.CNhat,
                                                  np.dot(int(np.floor_divide(x, self.tau)) *
                                                        (np.dot(self.Z.d2ANhat_d2q1,
                                                                np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x,self.tau))-1)) +
                                                                (np.floor_divide(x, self.tau) - 1) * np.dot(self.Z.dANhat_dq1,
                                                                                                            np.linalg.matrix_power(self.Z.ANhat,
                                                                                                                                   int(np.floor_divide(x, self.tau)) - 2)
                                                                                                            )
                                                        ),
                                                        self.Z.BNhat)
                                                 )

    # Convolution kernel
    def compute_L_i(self,theta,f,i):
        import numpy as np
        i+=1 # indexing in Python for i starts at 0 but we need it to start at 1. All signals have value 0 at t=0.
        if type(f) is np.ndarray:
            fsamp = f
        else:
            fsamp = np.array([f(x) for x in np.arange(1, i*self.tau,self.tau)])
        kernelsamp = np.array([self.conv_kernel(x) for x in np.arange(1, i * self.tau, self.tau)]).reshape(i)
        return np.sum([fsamp[j]*kernelsamp[i-j-2] for j in range(i)])

    def convolve_eifs(self,theta):
        convolved_eifs = np.zeros((self.P,self.n))
        for i in range(self.n):
            for j in range(self.P):
                e_j = np.zeros(self.P)
                e_j[j] = 1
                convolved_eifs[j, i] = self.compute_L_i(theta, self.g.eifs[j], i)
        return convolved_eifs

    # MCMC help
    def p_y_given_theta(self,theta,y):
        from scipy.stats import multivariate_normal
        v = 0.5*(self.vy_th + self.vy_th.T)
        p_y_given_th = multivariate_normal.pdf(self.y, np.zeros(self.n), v)
        return p_y_given_th

    def p_theta(self,theta):
        from scipy.stats import truncnorm
        return truncnorm.pdf(theta[1], 0, np.Infinity, 0.0046, 0.0006) * truncnorm.pdf(theta[2],0, np.Infinity, 1.23, 0.12)

    def acceptance(self,thetatry,thetaprev,y):
        p_thetatry = self.p_theta(thetatry)
        p_thetaprev = self.p_theta(thetaprev)
        num = self.p_y_given_theta(thetatry,y) * p_thetatry
        den = self.p_y_given_theta(thetaprev,y) * p_thetaprev
        return min(1, num / den)

    def EV_aP_given_theta_y(self,y):
        EV = [None, None]
        EV[0] = np.dot(np.dot(np.diag(self.g.eivs), self.lmatrix.T / self.vy_th), y.T)
        EV[1] = np.diag(self.g.eivs) - np.dot(np.dot(np.dot(np.diag(self.g.eivs), self.lmatrix.T), (np.linalg.solve(self.vy_th, self.lmatrix))), np.diag(self.g.eivs))
        EV[1] = (EV[1] + EV[1].T) / 2
        return EV

    # Set values of matrices needed for MCMC
    def set_d012_lmatrix(self,theta):
        self.lmatrix = self.convolve_eifs(theta)
        dkernsamp = np.array([self.dconv_kernel_th1(k*self.tau) for k in range(self.n)]).reshape(self.n)
        d2kernsamp = np.array([self.d2conv_kernel_th1(k*self.tau) for k in range(self.n)]).reshape(self.n)
        for j in range(self.P):
            eifsjsamp = [self.g.eifs[j](k*self.tau) for k in range(self.n)]
            for i in range(self.n):
                self.dlmatrix[j, i, 1] = sum([eifsjsamp[k]*dkernsamp[i-k-1] for k in range(i)])
                self.d2lmatrix[j, i, 1, 1] = sum([eifsjsamp[k]*d2kernsamp[i-k-1] for k in range(i)])
        self.dlmatrix[:, :, 2] = self.lmatrix / theta[2]

        self.d2lmatrix[:, :, 1, 2] = self.dlmatrix[:, :, 1] / theta[2]
        self.d2lmatrix[:, :, 2, 1] = self.d2lmatrix[:, :, 1, 2]

    def set_vy_thu(self,theta,u):
        for i in range(self.n):
            self.vy_thu[i, i] = (0.1 * self.compute_L_i(theta,u,i)) ** 2

    def set_dvy_thu(self,theta,u):
        pass

    def set_d2vy_thu(self,theta,u):
        pass

    def set_vy_th(self,theta):
        for i in range(self.n):
            for k in range(i):
                self.vy_th[i, k] = np.sum(np.multiply(np.multiply(self.g.eivs, self.lmatrix[:, i]), self.lmatrix[:, k]))
                if (i != k):
                    self.vy_th[k, i] = self.vy_th[i, k]
        self.vy_th += self.vy_thu

    def set_dvy_th(self,theta):
        for i in range(self.n):
            for k in range(self.n):
                for r in range(self.nTheta):
                    self.dvy_th[i, k, r] = np.sum(np.multiply(np.multiply(self.g.deivs[:, r],
                                                                          self.lmatrix[:, i]),
                                                              self.lmatrix[:, k]) +
                                                  np.multiply(self.g.eivs,(np.multiply(self.dlmatrix[:, i, r],
                                                                                       self.lmatrix[:, k]) +
                                                                           np.multiply(self.lmatrix[:, i], self.dlmatrix[:, k, r]))
                                                              )
                                                  )
        self.dvy_th += self.dvy_thu

    def set_d2vy_th(self,theta):
        for i in range(self.n):
            for k in range(self.n):
                for s in range(self.nTheta):
                    for r in range(self.nTheta):
                        self.d2vy_th[i, k, s, r] = np.sum(
                            np.multiply(np.multiply(self.g.d2eivs[:, s, r], self.lmatrix[:, i]), self.lmatrix[:, k]) +
                            np.multiply(self.g.deivs[:, r],
                                        np.multiply(self.dlmatrix[:, i, s], self.lmatrix[:, k]) +
                                        np.multiply(self.lmatrix[:, i],
                                                    self.dlmatrix[:, k, s])) +
                            np.multiply(self.g.deivs[:, s],
                                        np.multiply(self.dlmatrix[:, i, r], self.lmatrix[:, k]) +
                                        np.multiply(self.lmatrix[:, i], self.dlmatrix[:, k, r])) +
                            np.multiply(self.g.eivs,
                                        np.multiply(self.d2lmatrix[:, i, s, r], self.lmatrix[:, k]) +
                                        np.multiply(self.lmatrix[:, i],self.d2lmatrix[:, k, s, r]) +
                                        np.multiply(self.dlmatrix[:, i, s], self.dlmatrix[:, k, r]) +
                                        np.multiply(self.dlmatrix[:, i, r],self.dlmatrix[:, k, s])))
        self.d2vy_th += self.d2vy_thu

    def set_d2logpy_th(self,theta):
        from numpy import trace as tr
        from numpy.linalg import solve as d
        from numpy import dot as m
        s2 = self.d2vy_th
        s1 = self.dvy_th
        s0 = self.vy_th

        for s in range(self.nTheta):
            for r in range(self.nTheta):
                self.d2logpy_th[s, r] = -0.5 * (
                    tr(d(s0,(s2[:,:, s, r]))) -
                    tr(m(d(s0, s1[:, :, s]), d(s0,s1[:, :, r])) +
                                                            m(m(self.y,(m(d(s0,(m(s1[:, :, s], d(s0, s1[:, :, r])) -
                                                                   s2[:, :, s, r] +
                                                                   m(s1[:, :, r], d(s0, s1[:, :, s])))), np.linalg.inv(s0)))),self.y.T)))

    def set_vhat(self,theta):
        self.set_d2logpy_th(theta)
        self.vhat = - np.linalg.inv(self.d2logpy_th - self.infoTheta)
        self.vhat = 0.5 * (self.vhat + self.vhat.T)

    def calculate_new_operators(self,theta):
        self.g.set_theta(theta)
        self.set_vy_thu(theta,self.u)
        self.set_dvy_thu(theta,self.u)
        self.set_d2vy_thu(theta,self.u)
        self.set_d012_lmatrix(theta)
        self.set_vy_th(theta)
        self.set_dvy_th(theta)
        self.set_d2vy_th(theta)
        self.set_d2logpy_th(theta)
        # self.set_vhat(theta)

    """
    def plot_mcmc_results(self):
        fks = cell(1, K - burnin)

        # pillonetto_bell/src/mcmc.m:136
        for i in arange(1, K - burnin).reshape(-1):
            fks[i] = f_from_a_eifs(a[:, i].T, eifs)  # pillonetto_bell/src/mcmc.m:138

        fL_fU_fM = confidence_limits(fks)  # pillonetto_bell/src/mcmc.m:141
        fL = fL_fU_fM[1]
        # pillonetto_bell/src/mcmc.m:142
        fU = fL_fU_fM[2]
        # pillonetto_bell/src/mcmc.m:143
        fM = fL_fU_fM[3]
        # pillonetto_bell/src/mcmc.m:144
        plot(fM[t], 'b')
        hold('on')
        plot(sampled_u, 'ko')
        plot(fL[t], 'r--')
        plot(fU[t], 'r--')
        # Scatterplot thetas figure
        scatter(thetas[1, burnin:end()], thetas[2, burnin:end()], '.')
    """


if __name__ == "__main__":

    from scipy.stats import norm, truncnorm, multivariate_normal
    import numpy as np

    pb = PB()

    print "Calculating relevant operators...\n"
    pb.calculate_new_operators(pb.thetas[:, 0])
    pb.set_vhat(pb.thetas[:, 0])
    print "Finished calculating operators!\n"

    rejected = 0

    print "Beginning MCMC procedure\n"
    while pb.k < pb.K:

        if (pb.k % 500) == 0:
            print pb.k

        # Restricted to be nonnegative
        pb.thetas[:, pb.k] = pb.thetas[:, pb.k - 1] + np.sqrt(pb.alph) * np.dot(np.linalg.cholesky(pb.vhat), truncnorm.rvs(0,np.Infinity, size = pb.nTheta))

        for i in range(pb.n):
            pb.ys[i, pb.k] = pb.compute_L_i(pb.thetas[:, pb.k], pb.a[:, pb.k  - pb.burnin],i)
        pb.y = pb.ys[:, pb.k]

        try:
            acc=pb.acceptance(pb.thetas[:, pb.k], pb.thetas[:, pb.k - 1], pb.y)
        finally:
            pass

        c=np.random.uniform(0,1,1)
        if c > acc :
            if pb.k > pb.burnin:
                rejected += 1
            pb.thetas[:, pb.k]=pb.thetas[:, pb.k - 1]
        if pb.k < pb.burnin:
            pb.k += 1
            continue
        pb.recalculate_kernel(pb.thetas[:, pb.k])
        pb.calculate_new_operators(pb.thetas[:, pb.k])

        EV=pb.EV_aP_given_theta_y(pb.y)
        pb.a[:, pb.k + 1 - pb.burnin] = multivariate_normal(EV[0],EV[1])
        pb.k += 1

    # MCMC acceptance rate for theta
    print "Acceptance rate: " + 1 - rejected / (pb.k - pb.burnin)
