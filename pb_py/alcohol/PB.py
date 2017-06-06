class PB(object):

    def __init__(self):

        import numpy as np
        self.nTheta = 3
        self.P = 8
        self.burnin = 10
        self.K = 60
        self.n = 60
        self.tau = 0.1
        self.T = self.n * self.tau
        self.t = np.linspace(self.tau, self.T, self.n)
        self.alph = 3.

        # Actual sampled values. Attribute u only used for plotting.
        from scipy.stats import exponnorm
        self.y = np.array([exponnorm.pdf(x-3., 1.5) for x in np.arange(0.,self.T,self.tau)])
        self.u = np.concatenate(([1.],np.zeros(self.n-1)),axis=0)

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