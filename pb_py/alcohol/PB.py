class PB(object):

    def __init__(self):

        import numpy as np

        np.random.seed(990)

        # only used for plotting.
        self.u = np.array([0, 12.25, 24.5, 36.75, 49, 47, 45, 47.333, 49.667, 52, 49.25, 46.5, 43.75, 41, 39, 37, 36, 35, 34, 32.333, 30.667, 29, 27.667, 26.333, 25, 24.333, 23.667, 23, 22.167, 21.333, 20.5, 19.667, 18.833, 18, 17.333, 16.667, 16, 15, 14, 13, 12, 11, 10, 6.6667, 3.3333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])/100.
        self.y = np.array([0, 0, 0, 0, 1.7657, 6.1143, 11.289, 16.125, 19.03, 19.934, 20.982, 21.791, 24.744, 28.649, 32.077, 34.315, 34.315, 32.839, 32.553, 32.22, 31.887, 31.41, 31.077, 30.744, 30.601, 30.41, 30.458, 30.077, 29.887, 29.649, 29.22, 28.791, 28.553, 28.268, 28.125, 28.315, 28.172, 27.887, 27.649, 27.125, 26.934, 26.934, 26.649, 26.125, 25.41, 24.22, 23.22, 22.125, 21.268, 20.553, 19.791, 19.315, 18.506, 17.601, 16.744, 16.077, 15.315, 14.934, 14.553, 14.315, 13.839, 13.172, 12.363, 11.458, 10.934, 10.41, 10.03, 9.4105, 8.9343, 8.5533, 8.22, 7.8867, 7.5057, 6.9343, 6.4581, 6.2676, 6.1248, 6.22, 6.3152, 6.1724, 5.8867, 5.5533, 5.2676, 5.1248, 5.3152, 5.1724, 4.8867, 4.5533, 4.2676, 4.22, 4.1724, 3.8867, 3.5533, 3.2676, 3.22, 3.1724, 2.8867, 2.5533, 2.2676, 2.22, 2.1724, 1.8867, 1.5533, 1.2676, 1.1248, 1.22, 1.22, 1.22, 1.3152, 1.1724, 0.88667, 0.55333, 0.26762, 0.12476, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.24095, 0.20952, 0.14667, 0.073333, 0.010476, 0, 0, 0])/100.

        self.u = self.u[::2]
        self.y = self.y[::2]

        self.pad = 1

        self.nTheta = 3
        self.P = 6
        self.burnin = 4
        self.K = 12
        self.n = self.u.shape[0]
        self.tau = 5.
        self.T = self.n * self.tau
        self.t = np.linspace(self.tau, self.T, self.n)
        self.alph = 0.55

        # from scipy.stats import exponnorm
        # self.u = np.array([exponnorm.pdf((x - 3.) / 4., 0.3) for x in np.arange(self.tau, (self.n + 1) * self.tau, self.tau)])

        self.thetas = np.zeros((self.nTheta, self.K))
        self.thetas[:, 0] = [1., .0046, 1.23]

        self.a = np.zeros((self.P, self.K))
        self.ys = np.zeros((self.n, self.K))
        self.us = np.zeros((self.n, self.K))

        self.mats = None

    # MCMC help
    def p_y_given_theta(self, theta, y):
        from scipy.stats import multivariate_normal
        import numpy as np
        return multivariate_normal.pdf(y, mean=np.zeros(self.n), cov=self.mats.compute_vy_th(theta))

    def log_p_y_given_theta(self, theta, y):
        from scipy.stats import multivariate_normal
        import numpy as np
        return multivariate_normal.logpdf(y, mean=np.zeros(self.n), cov=self.mats.compute_vy_th(theta))

    def log_p_theta(self,theta):
        from scipy.stats import truncnorm
        return int(theta[0]>0 and theta[1]>0 and theta[2]>0) * truncnorm.logpdf(theta[1], a=0, b=10**10, loc=0.0046, scale=6e-4) * truncnorm.logpdf(theta[2], a=0, b=10**10, loc=1.23, scale=0.12)

    def acceptance(self,thetaprev,thetatry):
        from numpy import exp
        log_p_thetatry = self.log_p_theta(thetatry)
        log_p_thetaprev = self.log_p_theta(thetaprev)
        print "thetas:",thetaprev,thetatry
        print "log of priors:",log_p_thetaprev,log_p_thetatry
        den = self.log_p_y_given_theta(thetaprev, self.y)
        num = self.log_p_y_given_theta(thetatry, self.y)
        print "log( p(y|th) ):",den,num
        den += log_p_thetaprev
        num += log_p_thetatry
        return min(1, exp(num - den))

    def mean_cov_aP_given_theta_y(self,theta,y):
        import numpy as np
        self.mats.set_theta_update_operators(theta,self.u) # Here the u should not matter
        mean = np.dot(np.diag(self.mats.g.eivs), np.dot(self.mats.lmatrix, np.linalg.solve(self.mats.vy_th, y)))
        cov = np.diag(self.mats.g.eivs) - np.dot(np.dot(np.dot(np.diag(self.mats.g.eivs), self.mats.lmatrix),
                                                            np.linalg.solve(self.mats.vy_th, np.transpose(self.mats.lmatrix))),
                                                     np.diag(self.mats.g.eivs))
        cov = (cov + np.transpose(cov)) / 2
        return mean, cov

    def u_from_amp(self,amp):
        import numpy as np
        assert(amp.shape[0] == self.P)
        return np.dot(amp, np.array([self.mats.g.eifs[j](self.t) for j in range(self.P)]))

    def mcmc_init(self):
        from pbMatrices import pbMatrices
        import numpy as np
        print "Calculating matrices...\n"
        if self.y is None:
            self.mats = pbMatrices(self.thetas[:, 0], self.u, self.u, self.P, self.T, self.n, self.tau, k - 1)
            self.y = np.array([self.mats.compute_L_i(self.thetas[:, 0], self.u, i) for i in range(self.n)])
        else:
            self.mats = pbMatrices(self.thetas[:, 0], self.us[:,0], self.y, self.P, self.T, self.n, self.tau, k)
        self.vhat = self.mats.vhat
        print "alpha*V^ = ", self.alph*self.vhat
        print "Finished calculating matrices!\n"

    def mcmc(self):
        from scipy.stats import multivariate_normal
        import numpy as np
        import cPickle as pickle
        import gc

        self.mcmc_init()

        print "Beginning MCMC procedure\n"

        rejected = 0.
        # What step of the MCMC are we on?
        k = 1

        while k < self.K:

            print "\nWorking on step", k
            assert(all(x > 0. for x in self.thetas[:, k-1]))
            print np.linalg.eigvalsh(np.multiply(self.alph,self.vhat))
            assert(all(x>=0. for x in np.linalg.eigvalsh(np.multiply(self.alph,self.vhat))))

            # Restricted to be nonnegative
            self.thetas[:, k] = np.array(multivariate_normal.rvs(mean=self.thetas[:,k-1],
                                                                 cov=self.alph*self.vhat))

            print "this theta:", self.thetas[:,k]
            print "previous theta:", self.thetas[:,k-1]

            if any(x <= 0. for x in self.thetas[:,k]):
                acc = 0.
            else:
                acc = self.acceptance(self.thetas[:,k-1],self.thetas[:,k])

            c = np.random.uniform(0., 1.)
            if c > acc:
                self.thetas[:, k] = self.thetas[:,k-1]
                print "Rejected!", c, ">", acc
                print "theta reverted to", self.thetas[:, k-1]
                if k > self.burnin:
                    rejected += 1.

            mean, cov = self.mean_cov_aP_given_theta_y(self.thetas[:,k],self.y)

            self.a[:, k] = multivariate_normal.rvs(mean=mean, cov=cov)

            print "Amplitudes for draw", k, ":", self.a[:, k]
            self.us[:, k] = self.u_from_amp(self.a[:, k])

            for i in range(1,self.n+1):
                self.ys[i-1, k] = self.mats.compute_L_i(self.thetas[:,k], self.us[:, k], i)

            print "Completed step", k
            gc.collect()
            k += 1

        print "MCMC complete."
        # MCMC acceptance rate for theta
        print "Acceptance rate:", 1-rejected/(k-self.burnin)

        pickle.dump(self.us, open("us.pkl", "wb"))
        pickle.dump(self.ys, open("ys.pkl", "wb"))
        pickle.dump(self.thetas, open("thetas.pkl", "wb"))

    def plot_mcmc_results(self, loaded = True):
        import matplotlib.pyplot as plt
        import cPickle as pickle

        if loaded:
            self.us = pickle.load(open("us.pkl", "rb"))
            self.ys = pickle.load(open("ys.pkl", "rb"))
            self.thetas = pickle.load(open("thetas.pkl", "rb"))

        fL, fM, fU = self.confidence_limits()
        print fM
        fig1, ax1 = plt.subplots()
        if self.u is not None:
            ax1.plot(self.t, self.u, 'ko')
        if self.y is not None:
            ax1.plot(self.t, self.y, 'go')
        ax1.plot(self.t, fM, color='b')
        ax1.plot(self.t, fL, 'r--')
        ax1.plot(self.t, fU, 'r--')
        # Scatterplot thetas figure
        fig2, ax2 = plt.subplots()
        ax2.scatter(self.thetas[1, self.burnin:], self.thetas[2, self.burnin:])
        plt.show()

    def confidence_limits(self):
        from numpy import percentile, mean
        fL = percentile(self.us[:, self.burnin:], q=2, axis=1)
        fU = percentile(self.us[:, self.burnin:], q=98, axis=1)
        fM = mean(self.us[:,self.burnin:], axis=1)
        return fL, fM, fU

if __name__=="__main__":
    pb = PB()
    pb.mcmc()
    pb.plot_mcmc_results(True)
