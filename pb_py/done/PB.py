



class PB(object):

    def __init__(self):

        self.P = 20
        self.T = 1
        self.n = 50
        self.tau = 1. / 50.
        self.burnin = 400
        self.K = 8000
        self.alph = 3333
        self.nTheta = 2

        # What step of the MCMC are we on?
        self.k = 1

        import numpy as np
        self.infoTheta = np.array([ [0,0], [0,0.25] ])
        self.thetas = np.zeros((self.nTheta, self.K))
        self.thetas[:, 0] = [1, 10]
        self.t = np.linspace(self.tau, self.T, self.n)
        self.a = np.zeros((self.P, self.K - self.burnin))

        # Important matrices to be calculated
        self.vhat = np.zeros((self.nTheta,self.nTheta))
        self.lmatrix = np.zeros((self.P, self.n))
        self.dlmatrix = np.zeros((self.P,self.n,self.nTheta))
        self.d2lmatrix = np.zeros((self.P,self.n,self.nTheta,self.nTheta))
        self.vy_th = np.zeros((self.n,self.n))
        self.dvy_th = np.zeros((self.n, self.n, self.nTheta))
        self.d2vy_th = np.zeros((self.n,self.n,self.nTheta,self.nTheta))
        self.vy_thu = np.zeros((self.n,self.n))
        self.dvy_thu = np.zeros((self.n,self.n,self.nTheta))
        self.d2vy_thu = np.zeros((self.n,self.n,self.nTheta,self.nTheta))
        self.d2logpy_th = np.zeros(self.nTheta, self.nTheta)


        from rkhs import Green1_eigen
        self.g = Green1_eigen(self.P,self.T,self.thetas[:,0])

        self.u = None

        self.y = np.array([7.88231656866e-05, 0.00103141481843, 0.0042599302013, 0.0109414572813, 0.0216727450029, 0.0363139318338,
            0.0542174735416, 0.0742698619535, 0.0952090587957, 0.115750968435, 0.134715310806, 0.151123127707,
            0.164276892248, 0.173775284341, 0.179504813305, 0.181458101525, 0.180172458519, 0.175994829843,
            0.169610986984, 0.161628554675, 0.152658082896, 0.143436721543, 0.134590958132, 0.12662696759,
            0.119876351997, 0.11459774509, 0.110989979388, 0.108897237248, 0.10826521835, 0.10872970166, 0.10990919735,
            0.111354787636, 0.11254415309, 0.113050227885, 0.112295064502, 0.110246606493, 0.106490657257,
            0.101095550495, 0.094295603478, 0.0862019988128, 0.0773258488424, 0.0680937375501, 0.0589552679358,
            0.050283488449, 0.0424141229229, 0.0355148857564, 0.0296078491161, 0.0246289666774, 0.0204823531847,
            0.0170388615128])

    def set_u(self,u=None):
        if u is not None:
            self.u = u
        else:
            from scipy.stats import beta, norm
            from sympy import symbols
            s = symbols("s")
            self.u = norm.rvs(loc = 0.3*beta.pdf(s,12,7) + 0.6*beta.pdf(s,4,11), scale = 0.05)

    def set_y(self,y):
        self.y = y

    # Convolution kernel
    def compute_L_i(self,theta,f,i):
        from sympy import symbols, integrate
        from math import exp
        i+=1 # indexing in Python for i starts at 0 but we need it to start at 1
        x = symbols("x")
        h = exp(-theta[2] * (i*self.tau - x))
        return integrate(f*h, (x,0,i*self.tau))

    def convolve_eifs(self,theta):
        import numpy as np
        convolved_eifs = np.zeros((self.P,self.n))
        for i in range(self.n):
            for j in range(self.P):
                e_j = np.zeros(self.P)
                e_j[j] = 1
                convolved_eifs[j, i] = self.compute_L_i(theta, self.g.eifs[j], i)
        return convolved_eifs

    # MCMC help

    def p_y_given_theta(self,theta):
        from numpy.random import multivariate_normal
        self.set_vy_th(theta)
        p_y_given_th = multivariate_normal(self.y, 0, (self.vy_th + self.vy_th.T) / 2)
        return p_y_given_th

    def p_theta(self,theta):
        from scipy.stats import norm
        return int(theta[0] > 0 and theta[1] > 0) * norm.pdf(theta[1], 10, 2) / (1 - norm.cdf(- 5, 0, 1))

    def thetaprior(self,theta):
        from scipy.stats import norm
        return norm.pdf(theta, 0.0046, 0.0006)

    def acceptance(self,thetatry,thetaprev):
        p_thetatry = self.p_theta(thetatry)
        p_thetaprev = self.p_theta(thetaprev)
        num = self.p_y_given_theta(thetatry) * p_thetatry
        den = self.p_y_given_theta(thetaprev) * p_thetaprev
        return min(1, num / den)

    def EV_aP_given_theta_y(self):
        import numpy as np
        EV = list(2)
        EV[0] = np.dot(np.dot(np.diag(self.g.eivs), self.lmatrix.T / self.vy_th), self.y.T)
        EV[1] = np.diag(self.g.eivs) - np.dot(np.dot(np.dot(np.diag(self.g.eivs), self.lmatrix.T), (np.linalg.solve(self.vy_th, self.lmatrix))), np.diag(self.g.eivs))
        EV[1] = (EV[1] + EV[1].T) / 2
        return EV


    # Set values of matrices needed for MCMC

    def set_d012_lmatrix(self,theta):
        self.lmatrix = self.convolve_eifs(theta)

        self.dlmatrix[:, :, 0] = self.lmatrix / theta[0]
        self.dlmatrix[:, :, 1] = - 2 * theta[1] * self.lmatrix

        self.d2lmatrix[:, :, 0, 1] = - 2 * theta[1] / theta[0] * self.lmatrix
        self.d2lmatrix[:, :, 1, 0] = self.d2lmatrix[:, :, 1, 2]
        self.d2lmatrix[:, :, 1, 1] = ((4 * theta[1] ** 2) - 2) * self.lmatrix

    def set_vy_th(self,theta):
        import numpy as np
        for i in range(self.n):
            for k in range(i):
                self.vy_th[i, k] = sum(np.multiply(np.multiply(self.g.eivs, self.lmatrix[:, i]), self.lmatrix[:, k]))
                if i != k:
                    self.vy_th[k, i] = self.vy_th[i, k]
        self.vy_th += self.vy_thu

    def set_dvy_th(self,theta):
        import numpy as np
        self.g.set_theta(theta)
        self.set_d012_lmatrix(theta)
        for i in range(self.n):
            for k in range(self.n):
                for r in range(self.nTheta):
                    self.dvy_th[i, k, r] = sum(np.multiply(np.multiply(self.g.deivs[:, r],
                                                                       self.lmatrix[:, i]),
                                                           self.lmatrix[:, k]) +
                                               np.multiply(self.g.eivs,(np.multiply(self.dlmatrix[:, i, r],
                                                                                    self.lmatrix[:, k]) +
                                                                        np.multiply(self.lmatrix[:, i], self.dlmatrix[:, k, r]))
                                                           )
                                               )

        self.dvy_th += self.dvy_thu

    def set_d2vy_th(self,theta):
        import numpy as np
        self.g.set_theta(theta)
        self.set_d012_lmatrix(theta)
        self.d2vy_th = np.zeros((self.n, self.n, self.nTheta, self.nTheta))

        for i in range(self.n):
            for k in range(self.n):
                for s in range(self.nTheta):
                    for r in range(self.nTheta):
                        self.d2vy_th[i, k, s, r] = sum(
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

    def set_vy_thu(self,theta,u):
        for i in range(self.n):
            self.vy_thu[i, i] = (0.05 * self.compute_L_i(theta,u,i)) ** 2

    def set_dvy_thu(self,theta,u):
        pass

    def set_d2vy_thu(self,theta,u):
        pass

    def set_d2logpy_th(self,theta):
        import numpy as np

        assert(self.y.shape[0] == 1)

        self.set_d2vy_th(theta)
        s2 = self.d2vy_th
        self.set_dvy_th(theta)
        s1 = self.dvy_th
        self.set_vy_th(theta)
        s0 = self.vy_th

        self.d2logpy_th = np.zeros(self.nTheta, self.nTheta)

        for s in range(self.nTheta):
            for r in range(self.nTheta):
                self.d2logpy_th[s, r] = - 0.5 * (np.trace(np.linalg.solve(s0, (s2[:, :, s, r]))) -
                                            np.trace(np.dot(np.linalg.solve(s0, (s1[:, :, s] / s0)), s1[:, :, r])) +
                                            np.dot(np.dot(self.y,
                                                          np.linalg.solve(s0, (np.dot(s1[:, :, s] / s0, s1[:, :, r]) -
                                                                              s2[:, :, s, r] + np.dot(s1[:, :, r] / s0, s1[:, :, s]))) / s0), self.y.T))

    def set_vhat(self,theta = [1,10]):
        from numpy.linalg import inv
        self.set_d2logpy_th(theta)
        self.vhat = - inv(self.d2logpy_th - self.infoTheta)
        self.vhat = 0.5 * (self.vhat + self.vhat.T)





if __name__ == "__main__":

    from scipy.stats import norm,truncnorm
    import numpy as np

    pb = PB()

    ## Loop
    while pb.k < pb.K:

        if logical_not(pb.k % 500):
            k

        # Restricted to be nonnegative

        pb.thetas[:, pb.k] = pb.thetas[:, pb.k - 1] + truncnorm.rvs(pb.alph * pb.vhat) # pillonetto_bell/src/mcmc.m:88
        try:
            acc=acceptance(thetas[:,k],thetas[:,k - 1] , y,tau,T,P,n, rkhs_eigenfile, data_path) # pillonetto_bell/src/mcmc.m:91
        finally:
            pass
        c=unifrnd(0,1,1) # pillonetto_bell/src/mcmc.m:100
        if c > acc :
            if k > burnin:
                rejected=rejected + 1
            # pillonetto_bell/src/mcmc.m:104
            thetas[:,k]=thetas[:,k - 1]
        # pillonetto_bell/src/mcmc.m:106
        if k < burnin:
            k=k + 1
            # pillonetto_bell/src/mcmc.m:110
            continue
        EV=EV_aP_given_theta_y(thetas[:,k],y,rkhs_eigenfile,P,T,n,tau, data_path) # pillonetto_bell/src/mcmc.m:114
        mus= EV [1 ]
        # pillonetto_bell/src/mcmc.m:116
        covs=EV[2] # pillonetto_bell/src/mcmc.m:117
        a[:,k + 1 - burnin]=mvnrnd( mus,covs)
        # pillonetto_bell/src/mcmc.m:119 #     progress.step = k; #     progress.theta1=[thetas(1,k-1),thetas(1,k)];
        #     progress.theta2=[thetas(2,k-1),thetas(2,k)];
        #     progress.acc1=acc(1);
        #     progress.acc2=acc(2);
        #     progress
        k=k + 1
    # pillonetto_bell/src/mcmc.m:129


    # MCMC acceptance rate for theta
    1 - rejected / (k - burnin)
    fks=cell(1,K - burnin)
    # pillonetto_bell/src/mcmc.m:136
    for i in arange(1,K - burnin).reshape(-1):
        fks[i]=f_from_a_eifs(a[:,i].T,eifs) # pillonetto_bell/src/mcmc.m:138

    fL_fU_fM=confidence_limits(fks) # pillonetto_bell/src/mcmc.m:141
    fL=fL_fU_fM[1]
    # pillonetto_bell/src/mcmc.m:142
    fU=fL_fU_fM[2]
    # pillonetto_bell/src/mcmc.m:143
    fM= fL_fU_fM[3]
    # pillonetto_bell/src/mcmc.m:144
    plot(fM [ t],'b')
    hold('on')
    plot(sampled_u,'ko') plot(fL[t ],'r--')
    plot(fU[t],'r--')
    # Scatterplot thetas figure
    scatter(thetas[1,burnin:end()],thetas [2,burnin:end()],'.')