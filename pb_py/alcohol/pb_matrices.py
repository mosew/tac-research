class pbMatrices(object):

    def __init__(self, theta, P, T, n, tau, k):
        self.theta = theta
        self.nTheta = theta.shape([0])
        self.P = P
        self.T = T
        self.n = n
        self.tau = tau
        self.k = k

        from Parabolic_System import Parabolic_System
        self.Z = Parabolic_System(self.theta[1],self.theta[2])
        self.Z.define_operators()

        import numpy as np
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

        self.conv_kernel = lambda x: np.dot(self.Z.CNhat, np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau))), self.Z.BNhat))
        self.dconv_kernel_th1 = lambda x: np.dot(self.Z.CNhat,
                                                 int(np.floor_divide(x, self.tau)) *
                                                 np.dot(self.Z.dANhat_dq1,
                                                        np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau)-1)),
                                                               self.Z.BNhat)))
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

        from rkhs import Green1_eigen
        self.g = Green1_eigen(self.P,self.T,self.theta)

    def compute_kernel(self,theta):
        from Parabolic_System import Parabolic_System
        import numpy as np

        conv_kernel = lambda x: np.dot(self.Z.CNhat, np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau))), self.Z.BNhat))
        dconv_kernel_th1 = lambda x: np.dot(self.Z.CNhat, int(np.floor_divide(x, self.tau)) * np.dot(self.Z.dANhat_dq1, np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau)-1)), self.Z.BNhat)))
        d2conv_kernel_th1 = lambda x: np.dot(self.Z.CNhat,
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
        return conv_kernel, dconv_kernel_th1, d2conv_kernel_th1

    def reset_kernel(self,theta):
        self.conv_kernel, self.dconv_kernel_th1, self.d2conv_kernel_th1 = self.compute_kernel(theta)

    # Convolution kernel
    def compute_L_i(self, theta, f, i):
        import numpy as np
        i+=1 # indexing in Python for i starts at 0 but we need it to start at 1. All signals have value 0 at t=0.
        if type(f) is np.ndarray:
            fsamp = f[:i]
        else:
            fsamp = np.array([f(x) for x in np.arange(1, i * self.tau, self.tau)])
        kernelsamp = np.array([self.conv_kernel(x) for x in np.arange(1, i * self.tau, self.tau)]).reshape(i)
        return np.sum([fsamp[j]*kernelsamp[i-j-1] for j in range(i)])

    def convolve_eifs(self,theta):
        import numpy as np
        convolved_eifs = np.zeros((self.P, self.n))
        for i in range(self.n):
            for j in range(self.P):
                convolved_eifs[j, i] = self.compute_L_i(theta, self.g.eifs[j], i)
        return convolved_eifs

    def f_from_amp(self,amp):
        import numpy as np
        return np.dot(amp, np.array([self.g.eifs[j](self.t) for j in range(amp.shape[0])]))


    # Compute / set values of matrices needed for MCMC
    def isthetaprev(self, theta):
        return all(x == 0 for x in theta - self.thetas[:, self.k-1])

    def compute_d012_lmatrix(self,theta):
        import numpy as np
        lmatrix = self.convolve_eifs(theta)
        dlmatrix = np.zeros(self.dlmatrix.shape)
        d2lmatrix = np.zeros(self.d2lmatrix.shape)
        dkernsamp = np.array([self.dconv_kernel_th1(k*self.tau) for k in range(self.n)]).reshape(self.n)
        d2kernsamp = np.array([self.d2conv_kernel_th1(k*self.tau) for k in range(self.n)]).reshape(self.n)
        for j in range(self.P):
            eifsjsamp = [self.g.eifs[j](k*self.tau) for k in range(self.n)]
            for i in range(self.n):
                dlmatrix[j, i, 1] = sum([eifsjsamp[k]*dkernsamp[i-k-1] for k in range(i)])
                d2lmatrix[j, i, 1, 1] = sum([eifsjsamp[k]*d2kernsamp[i-k-1] for k in range(i)])
        dlmatrix[:, :, 2] = lmatrix / theta[2]

        d2lmatrix[:, :, 1, 2] = dlmatrix[:, :, 1] / theta[2]
        d2lmatrix[:, :, 2, 1] = d2lmatrix[:, :, 1, 2]
        return lmatrix, dlmatrix, d2lmatrix

    def set_d012_lmatrix(self, theta):
        self.lmatrix, self.dlmatrix, self.d2lmatrix = self.compute_d012_lmatrix(theta)

    def compute_vy_thu(self,theta,u):
        import numpy as np
        vy_thu = np.zeros(self.vy_thu.shape)
        for i in range(self.n):
            vy_thu[i, i] = theta[0]**2
        return np.array(vy_thu)

    def set_vy_thu(self,theta,u):
        self.vy_thu = self.compute_vy_thu(theta,u)

    def compute_dvy_thu(self,theta,u):
        import numpy as np
        dvy_thu = np.zeros(self.dvy_thu.shape)
        for i in range(self.n):
            dvy_thu[i, i, 0] = 2*theta[0]
        return np.array(dvy_thu)

    def set_dvy_thu(self,theta,u):
        self.dvy_thu = self.compute_dvy_thu(theta,u)

    def compute_d2vy_thu(self,theta,u):
        import numpy as np
        d2vy_thu = np.zeros(self.d2vy_thu.shape)
        d2vy_thu[:, :, 0, 0] = 2 * np.identity(self.n)
        return d2vy_thu

    def set_d2vy_thu(self,theta,u):
        self.d2vy_thu = self.compute_d2vy_thu(theta,u)

    def compute_vy_th(self,theta):
        import numpy as np
        from rkhs import Green1_eigen

        vy_th = np.zeros(self.vy_th.shape)

        # This part could be made faster by exploiting simplicity of Green1 eigenvalues dependence on theta
        if any(x != 0 for x in (theta - self.g.theta)):
            g = Green1_eigen(self.P, self.T, theta)
        else:
            g = self.g

        if not self.isthetaprev(theta):
            lmatrix, _, _ = self.compute_d012_lmatrix(theta)
            vy_thu = self.compute_vy_thu(theta, self.us[:, 0]) # By hypothesis, vy_th can't depend on u...?
        else:
            lmatrix = self.lmatrix
            vy_thu = self.vy_thu

        for i in range(self.n):
            for k in range(i+1):
                vy_th[i, k] = np.sum(np.multiply(np.multiply(g.eivs, lmatrix[:, i]), lmatrix[:, k]))
                if i != k:
                    vy_th[k, i] = vy_th[i, k]
        vy_th += vy_thu
        return vy_th

    def set_vy_th(self,theta):
        self.vy_th = self.compute_vy_th(theta)

    def compute_dvy_th(self,theta):
        import numpy as np
        from rkhs import Green1_eigen
        dvy_th = np.zeros(self.dvy_th.shape)

        if any(x != 0 for x in (theta - self.g.theta)):
            g = Green1_eigen(self.P, self.T, theta)
        else:
            g = self.g

        if not self.isthetaprev(theta):
            lmatrix, dlmatrix, _ = self.compute_d012_lmatrix(theta)
            dvy_thu = self.compute_dvy_thu(theta, self.us[:, 0])
        else:
            lmatrix = self.lmatrix
            dlmatrix = self.dlmatrix
            dvy_thu = self.dvy_thu

        for i in range(self.n):
            for k in range(self.n):
                for r in range(self.nTheta):
                    dvy_th[i, k, r] = np.sum(np.multiply(np.multiply(g.deivs[:, r], lmatrix[:, i]),
                                                         lmatrix[:, k]) +
                                             np.multiply(g.eivs,(np.multiply(dlmatrix[:, i, r],
                                                                             lmatrix[:, k]) +
                                                                 np.multiply(lmatrix[:, i], dlmatrix[:, k, r]))
                                                         )
                                             )
        dvy_th += dvy_thu
        return dvy_th

    def set_dvy_th(self,theta):
        self.dvy_th = self.compute_dvy_th(theta)

    def compute_d2vy_th(self,theta):
        import numpy as np
        from rkhs import Green1_eigen

        d2vy_th = np.zeros(self.d2vy_th.shape)

        if any(x != 0 for x in (theta - self.g.theta)):
            g = Green1_eigen(self.P, self.T, theta)
        else:
            g = self.g

        if not self.isthetaprev(theta):
            lmatrix, dlmatrix, d2lmatrix = self.compute_d012_lmatrix(theta)
            d2vy_thu = self.compute_d2vy_thu(theta, self.us[:, 0])
        else:
            lmatrix = self.lmatrix
            dlmatrix = self.dlmatrix
            d2lmatrix = self.d2lmatrix
            d2vy_thu = self.d2vy_thu

        for i in range(self.n):
            for k in range(self.n):
                for s in range(self.nTheta):
                    for r in range(self.nTheta):
                        self.d2vy_th[i, k, s, r] = np.sum(
                            np.multiply(np.multiply(g.d2eivs[:, s, r], lmatrix[:, i]), lmatrix[:, k]) +
                            np.multiply(g.deivs[:, r],
                                        np.multiply(dlmatrix[:, i, s], lmatrix[:, k]) +
                                        np.multiply(lmatrix[:, i],
                                                    dlmatrix[:, k, s])) +
                            np.multiply(g.deivs[:, s],
                                        np.multiply(dlmatrix[:, i, r], lmatrix[:, k]) +
                                        np.multiply(lmatrix[:, i], dlmatrix[:, k, r])) +
                            np.multiply(g.eivs,
                                        np.multiply(d2lmatrix[:, i, s, r], lmatrix[:, k]) +
                                        np.multiply(lmatrix[:, i],d2lmatrix[:, k, s, r]) +
                                        np.multiply(dlmatrix[:, i, s], dlmatrix[:, k, r]) +
                                        np.multiply(dlmatrix[:, i, r], dlmatrix[:, k, s])))
        d2vy_th += d2vy_thu
        return d2vy_th

    def set_d2vy_th(self,theta):
        self.d2vy_th = self.compute_d2vy_th(theta)

    def set_d2logpy_th(self,theta):
        import numpy as np
        from numpy import trace as tr
        from numpy.linalg import solve as d
        from numpy import dot as m

        if not self.isthetaprev(theta):
            s2 = self.compute_d2vy_th(theta)
            s1 = self.compute_dvy_th(theta)
            s0 = self.compute_vy_th(theta)
        else:
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
        import numpy as np
        self.vhat = np.linalg.inv(- self.d2logpy_th + self.infoTheta)
        self.vhat = 0.5 * (self.vhat + self.vhat.T)

    def set_theta_update_operators(self,theta,u):
        self.g.set_theta(theta)
        self.set_vy_thu(theta,u)
        self.set_dvy_thu(theta,u)
        self.set_d2vy_thu(theta,u)
        self.set_d012_lmatrix(theta)
        self.set_vy_th(theta)
        self.set_dvy_th(theta)
        self.set_d2vy_th(theta)
        self.set_d2logpy_th(theta)
        self.set_vhat(theta)
