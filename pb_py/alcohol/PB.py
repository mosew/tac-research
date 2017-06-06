class PB(object):

    def __init__(self):

        import numpy as np
        # only used for plotting.
        self.u = np.array([0, 12.25, 24.5, 36.75, 49, 47, 45, 47.333, 49.667, 52, 49.25, 46.5, 43.75, 41, 39, 37, 36, 35, 34, 32.333, 30.667, 29, 27.667, 26.333, 25, 24.333, 23.667, 23, 22.167, 21.333, 20.5, 19.667, 18.833, 18, 17.333, 16.667, 16, 15, 14, 13, 12, 11, 10, 6.6667, 3.3333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.y = np.array([1.1895e-10, 0, 0, 0, 1.7657, 6.1143, 11.289, 16.125, 19.03, 19.934, 20.982, 21.791, 24.744, 28.649, 32.077, 34.315, 34.315, 32.839, 32.553, 32.22, 31.887, 31.41, 31.077, 30.744, 30.601, 30.41, 30.458, 30.077, 29.887, 29.649, 29.22, 28.791, 28.553, 28.268, 28.125, 28.315, 28.172, 27.887, 27.649, 27.125, 26.934, 26.934, 26.649, 26.125, 25.41, 24.22, 23.22, 22.125, 21.268, 20.553, 19.791, 19.315, 18.506, 17.601, 16.744, 16.077, 15.315, 14.934, 14.553, 14.315, 13.839, 13.172, 12.363, 11.458, 10.934, 10.41, 10.03, 9.4105, 8.9343, 8.5533, 8.22, 7.8867, 7.5057, 6.9343, 6.4581, 6.2676, 6.1248, 6.22, 6.3152, 6.1724, 5.8867, 5.5533, 5.2676, 5.1248, 5.3152, 5.1724, 4.8867, 4.5533, 4.2676, 4.22, 4.1724, 3.8867, 3.5533, 3.2676, 3.22, 3.1724, 2.8867, 2.5533, 2.2676, 2.22, 2.1724, 1.8867, 1.5533, 1.2676, 1.1248, 1.22, 1.22, 1.22, 1.3152, 1.1724, 0.88667, 0.55333, 0.26762, 0.12476, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.24095, 0.20952, 0.14667, 0.073333, 0.010476, 0, 1.6217e-10, 2.5312e-12])
        self.pad = 8

        self.u = np.concatenate((np.zeros(self.pad), self.u))
        self.y = np.concatenate((np.zeros(self.pad), self.y))

        self.nTheta = 3
        self.P = 8
        self.burnin = 10
        self.K = 60
        self.n = self.u.shape[0]
        self.tau = 1./12.
        self.T = self.n * self.tau
        self.t = np.linspace(self.tau, self.T, self.n)
        self.alph = 1000000.

        # What step of the MCMC are we on?
        self.k = 1

        self.thetas = np.zeros((self.nTheta, self.K))
        self.thetas[:, 0] = [1., .0046, 1.23]

        self.a = np.zeros((self.P, self.K))
        self.ys = np.zeros((self.n, self.K))
        self.us = np.zeros((self.n, self.K))

        self.mats_try = None
        self.mats_prev = None

    def mcmc_init(self):
        from pbMatrices import pbMatrices
        print "Calculating matrices...\n"
        self.mats_prev = pbMatrices(self.thetas[:, 0], self.u, self.ys[:, 0], self.P, self.T, self.n, self.tau, self.k - 1)
        self.y = np.array([self.mats_prev.compute_L_i(self.thetas[:, 0], self.u, i) for i in range(self.n)])
        self.mats_prev = pbMatrices(self.thetas[:, 0], self.u, self.y, self.P, self.T, self.n, self.tau, self.k - 1)
        self.mats_try = self.mats_prev
        print "Finished calculating matrices!\n"

    def mcmc_draw_and_accept(self):
        pass

    # MCMC help
    def p_y_given_theta(self, theta, y, prev_flag):
        from scipy.stats import multivariate_normal
        import numpy as np
        if prev_flag:
            return multivariate_normal.pdf(y, mean=np.zeros(self.n), cov=self.mats_prev.vy_th)
        else:
            return multivariate_normal.pdf(y, mean=np.zeros(self.n), cov=self.mats_try.vy_th)

    def p_theta(self,theta):
        from scipy.stats import truncnorm
        import numpy as np
        return truncnorm.pdf(theta[1], 0, np.Infinity, 0.0046, 0.0006) * truncnorm.pdf(theta[2],0, np.Infinity, 1.23, 0.12)

    def acceptance(self,y):
        thetaprev = self.mats_prev.theta
        thetatry = self.mats_try.theta
        p_thetatry = self.p_theta(thetatry)
        p_thetaprev = self.p_theta(thetaprev)
        den = self.p_y_given_theta(thetaprev, y, True) * p_thetaprev
        num = self.p_y_given_theta(thetatry, y, False) * p_thetatry
        return min(1, num / den)

    def mean_cov_aP_given_theta_y(self,y):
        import numpy as np
        mean = np.dot(np.diag(self.mats_try.g.eivs), np.dot(self.mats_try.lmatrix, np.linalg.solve(self.mats_try.vy_th, y)))
        cov = np.diag(self.mats_try.g.eivs) - np.dot(np.dot(np.dot(np.diag(self.mats_try.g.eivs), self.mats_try.lmatrix),
                                                            np.linalg.solve(self.mats_try.vy_th, np.transpose(self.mats_try.lmatrix))),
                                                     np.diag(self.mats_try.g.eivs))
        cov = (cov + np.transpose(cov)) / 2
        return mean, cov

    def u_from_amp(self,amp):
        import numpy as np
        return np.dot(amp, np.array([self.mats_try.g.eifs[j](self.t) for j in range(amp.shape[0])]))

    def mcmc(self):
        from scipy.stats import multivariate_normal
        from pbMatrices import pbMatrices
        import numpy as np
        import cPickle as pickle

        self.mcmc_init()

        print "Beginning MCMC procedure\n"

        rejected = 0.

        while self.k < self.K:
            print "Working on step ", self.k
            assert(all(x > 0. for x in self.thetas[:, self.k-1]))

            # Restricted to be nonnegative
            self.thetas[:, self.k] = np.array(multivariate_normal.rvs(mean=self.mats_prev.theta,
                                                                      cov=np.sqrt(self.alph)*self.mats_prev.vhat))

            self.mats_try = pbMatrices(self.thetas[:, self.k], self.us[:, self.k - 1], self.y, self.P, self.T, self.n, self.tau, self.k)

            if any(x <= 0. for x in self.mats_try.theta):
                acc = 0.
            else:
                print "Calculating acceptance ratio."
                acc=self.acceptance(self.y)

            c = np.random.uniform(0., 1., 1)
            if c > acc:
                self.thetas[:, self.k] = self.mats_prev.theta
                self.mats_try = self.mats_prev
                if self.k > self.burnin:
                    rejected += 1.

            if self.k < self.burnin:
                self.k += 1
                continue

            mean, cov = self.mean_cov_aP_given_theta_y(self.y)

            self.a[:, self.k] = multivariate_normal.rvs(mean=mean, cov=cov)
            self.us[:, self.k] = self.u_from_amp(self.a[:, self.k])

            for i in range(self.n):
                self.ys[i, self.k] = self.mats_try.compute_L_i(self.mats_try.theta, self.us[:, self.k], i)

            self.k += 1
            self.mats_prev = self.mats_try
            print "Completed step k = ", self.k

        print "MCMC complete."
        # MCMC acceptance rate for theta
        print "Acceptance rate: ", 1-rejected/(self.k-self.burnin)

        pickle.dump(self.us, open("us.pkl", "wb"))
        pickle.dump(self.ys, open("ys.pkl", "wb"))
        pickle.dump(self.thetas, open("thetas.pkl", "wb"))

    def plot_mcmc_results(self, loaded = False):
        import matplotlib.pyplot as plt
        import cPickle as pickle

        if loaded:
            self.us = pickle.load(open("us.pkl", "rb"))
            self.ys = pickle.load(open("ys.pkl", "rb"))
            self.thetas = pickle.load(open("thetas.pkl", "rb"))

        fL, fM, fU = self.confidence_limits()
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
    from scipy.stats import exponnorm
    import numpy as np

    pb = PB()
    pb.mcmc()
    pb.plot_mcmc_results(True)