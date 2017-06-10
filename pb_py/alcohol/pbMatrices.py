class pbMatrices(object):

    def __init__(self, theta, u, y, P, T, n, tau, k, vhat=None):
        self.theta = theta
        self.nTheta = theta.shape[0]
        self.P = P
        self.T = T
        self.n = n
        self.tau = tau
        self.k = k
        self.y = y
        self.u = u

        from Parabolic_System import Parabolic_System
        self.Z = Parabolic_System(self.theta[1],self.theta[2], self.u, self.n, self.tau)
        self.Z.define_operators()

        import numpy as np
        self.infoTheta = np.array([[0.,0.,0.], [0.,(.0006)**(-2), 0.],[0., 0.,(0.12)**(-2)]])
        # Important matrices to be calculated
        if vhat is not None:
            self.vhat = vhat
        else:
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
        self.g = Green1_eigen(self.P,self.T,self.theta)
        self.set_theta_update_operators(self.theta, self.u)

    def conv_kernel(self, x):
        import numpy as np
        assert(x==float(x))
        if int(x / self.tau)==0:
            return 0.
        else:
            return np.dot(self.Z.CNhat, np.dot(np.linalg.matrix_power(self.Z.ANhat,
                                                                      int(x/self.tau)-1),
                                               self.Z.BNhat))

    def dconv_kernel_th1(self, x):
        import numpy as np
        if int(np.floor_divide(x, self.tau)) <= 1:
            return 0.
        else:
            return np.dot(self.Z.CNhat,
                          np.multiply(int(np.floor_divide(x, self.tau)-1.),
                          np.dot(self.Z.dANhat_dq1,
                                 np.dot(np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x, self.tau)-2.)),
                                        self.Z.BNhat))))

    def d2conv_kernel_th1(self, x):
        import numpy as np
        if int(np.floor_divide(x, self.tau)) <= 2:
            return 0.
        else:
            return np.dot(self.Z.CNhat,
                          np.dot(int(np.floor_divide(x, self.tau)-1) *
                                 (np.dot(self.Z.d2ANhat_d2q1,
                                         np.linalg.matrix_power(self.Z.ANhat, int(np.floor_divide(x,self.tau))-2)) +
                                  (np.floor_divide(x, self.tau) - 2) * np.dot(self.Z.dANhat_dq1,
                                                                              np.linalg.matrix_power(self.Z.ANhat,
                                                                                                     int(np.floor_divide(x, self.tau)) - 3)
                                                                              )
                                  ),
                                 self.Z.BNhat)
                          )

    def reset_kernel(self,theta):
        self.Z.reset_q_u(theta[1],theta[2],self.u)

    # Convolution kernel
    def compute_L_i(self, theta, f, i):
        import numpy as np
        i+=1 # indexing in Python for i starts at 0 but we need it to start at 1. All signals have value 0 at t=0.
        ti = np.linspace(self.tau, i*self.tau, i)
        if type(f) is np.ndarray:
            fsamp = np.array(f[:i])
        else:
            fsamp = np.array([f(x) for x in ti])
        kernelsamp = np.array([self.conv_kernel(x) for x in ti]).reshape(i)
        return np.sum([fsamp[j]*kernelsamp[i-j-1] for j in range(i)])

    def convolve_eifs(self,theta):
        import numpy as np
        convolved_eifs = np.zeros((self.P, self.n))
        for i in range(self.n):
            for j in range(self.P):
                convolved_eifs[j, i] = self.compute_L_i(theta, self.g.eifs[j], i)
        return convolved_eifs

    # Compute / set values of matrices needed for MCMC
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

        d2lmatrix[:, :, 1, 2] = dlmatrix[:, :, 2] / theta[2]
        d2lmatrix[:, :, 2, 1] = d2lmatrix[:, :, 1, 2]
        return lmatrix, dlmatrix, d2lmatrix

    def set_d012_lmatrix(self, theta):
        self.lmatrix, self.dlmatrix, self.d2lmatrix = self.compute_d012_lmatrix(theta)

    def compute_vy_thu(self,theta):
        import numpy as np
        vy_thu = np.zeros(self.vy_thu.shape)
        for i in range(self.n):
            vy_thu[i, i] = (theta[0]*theta[2])**2
        return np.array(vy_thu)

    def set_vy_thu(self,theta):
        self.vy_thu = self.compute_vy_thu(theta)

    def compute_dvy_thu(self,theta):
        import numpy as np
        dvy_thu = np.zeros(self.dvy_thu.shape)
        for i in range(self.n):
            dvy_thu[i, i, 0] = 2*theta[0]*(theta[2]**2)
            dvy_thu[i, i, 2] = 2*theta[2]*(theta[0]**2)
        return np.array(dvy_thu)

    def set_dvy_thu(self,theta):
        self.dvy_thu = self.compute_dvy_thu(theta)

    def compute_d2vy_thu(self,theta):
        import numpy as np
        d2vy_thu = np.zeros(self.d2vy_thu.shape)
        d2vy_thu[:, :, 0, 0] = np.multiply(2*(theta[2]**2), np.identity(self.n))
        d2vy_thu[:, :, 2, 2] = np.multiply(2*(theta[0]**2), np.identity(self.n))
        d2vy_thu[:, :, 0, 2] = np.multiply(4*theta[0]*theta[2], np.identity(self.n))
        d2vy_thu[:, :, 2, 0] = d2vy_thu[:, :, 0, 2]
        return d2vy_thu

    def set_d2vy_thu(self,theta):
        self.d2vy_thu = self.compute_d2vy_thu(theta)

    def compute_vy_th(self,theta):
        import numpy as np
        vy_th = np.zeros(self.vy_th.shape)
        for i in range(self.n):
            for k in range(i+1):
                vy_th[i, k] = np.sum(np.multiply(np.multiply(self.g.eivs, self.lmatrix[:, i]), self.lmatrix[:, k]))
                if i != k:
                    vy_th[k, i] = vy_th[i, k]
        vy_th += self.vy_thu
        vy_th = 0.5*(vy_th + np.transpose(vy_th))
        return vy_th

    def set_vy_th(self,theta):
        self.vy_th = self.compute_vy_th(theta)

    def compute_dvy_th(self,theta):
        import numpy as np
        dvy_th = np.zeros(self.dvy_th.shape)

        for i in range(self.n):
            for k in range(self.n):
                for r in range(self.nTheta):
                    dvy_th[i, k, r] = np.sum(np.multiply(np.multiply(self.g.deivs[:, r], self.lmatrix[:, i]),
                                                         self.lmatrix[:, k]) +
                                             np.multiply(self.g.eivs,(np.multiply(self.dlmatrix[:, i, r],
                                                                                  self.lmatrix[:, k]) +
                                                                      np.multiply(self.lmatrix[:, i], self.dlmatrix[:, k, r]))
                                                         )
                                             )
        dvy_th += self.dvy_thu
        return dvy_th

    def set_dvy_th(self,theta):
        self.dvy_th = self.compute_dvy_th(theta)

    def compute_d2vy_th(self,theta):
        import numpy as np

        d2vy_th = np.zeros(self.d2vy_th.shape)

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
                                        np.multiply(self.lmatrix[:, i], self.d2lmatrix[:, k, s, r]) +
                                        np.multiply(self.dlmatrix[:, i, s], self.dlmatrix[:, k, r]) +
                                        np.multiply(self.dlmatrix[:, i, r], self.dlmatrix[:, k, s])))
        d2vy_th += self.d2vy_thu
        return d2vy_th

    def set_d2vy_th(self,theta):
        self.d2vy_th = self.compute_d2vy_th(theta)

    def set_d2logpy_th(self,theta):
        import numpy as np
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
                                      m(s1[:, :, r], d(s0, s1[:, :, s])))), np.linalg.inv(s0)))),np.transpose(self.y))))

    def set_vhat(self,theta):
        import numpy as np
        self.vhat = np.linalg.inv(self.d2logpy_th - self.infoTheta)
        self.vhat = 0.5 * (self.vhat + np.transpose(self.vhat))

    def set_theta_update_operators(self,theta, u):
        self.theta = theta
        self.Z.reset_q_u(theta[1],theta[2],u)
        self.g.set_theta(theta)
        self.set_vy_thu(theta)
        self.set_dvy_thu(theta)
        self.set_d2vy_thu(theta)
        self.set_d012_lmatrix(theta)
        self.set_vy_th(theta)
        self.set_dvy_th(theta)
        self.set_d2vy_th(theta)
        self.set_d2logpy_th(theta)
        self.set_vhat(theta)
