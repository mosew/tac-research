
class PB(object):

    def __init__(self,parabolic_system):
        self.data_path = None
        self.P = 20
        self.T = 1
        self.n = 50
        self.tau = 1 / 50
        tau = parabolic_system.tau
        self.t = range(tau, self.T, tau)
        self.burnin = 400
        self.K = 8000
        self.alph = 3333
        self.nTheta = 2

        import numpy as np
        self.thetas = np.zeros((self.nTheta, self.K))
        self.thetas[:, 0] = [1, 10]
        self.a = np.zeros((self.P, self.K - self.burnin))

        from rkhs import Green1_eigen
        g = Green1_eigen(self.P,self.T,self.thetas[:,0])

        self.y = np.array([7.88231656866e-05, 0.00103141481843, 0.0042599302013, 0.0109414572813, 0.0216727450029, 0.0363139318338,
            0.0542174735416, 0.0742698619535, 0.0952090587957, 0.115750968435, 0.134715310806, 0.151123127707,
            0.164276892248, 0.173775284341, 0.179504813305, 0.181458101525, 0.180172458519, 0.175994829843,
            0.169610986984, 0.161628554675, 0.152658082896, 0.143436721543, 0.134590958132, 0.12662696759,
            0.119876351997, 0.11459774509, 0.110989979388, 0.108897237248, 0.10826521835, 0.10872970166, 0.10990919735,
            0.111354787636, 0.11254415309, 0.113050227885, 0.112295064502, 0.110246606493, 0.106490657257,
            0.101095550495, 0.094295603478, 0.0862019988128, 0.0773258488424, 0.0680937375501, 0.0589552679358,
            0.050283488449, 0.0424141229229, 0.0355148857564, 0.0296078491161, 0.0246289666774, 0.0204823531847,
            0.0170388615128])

    def set_y(self,y):
        self.y = y

    def plot_u_vs_y(self,u,y):
        from matplotlib import pyplot as plt


    # plot(t,y)

    ## Generate parameter variance estimates
    # I think technically this should be re-evaluated at the current draw of theta,
    # but I don't think it should vary all that much.
    Vhat_=Vhat( y ,cat(1, 10).T, tau,P, T,n, eivs, rkhs_eigenfile,data_path )
    # pillonetto_bell/src/mcmc.m:72
    Vhat_=(Vhat_ + Vhat_.T) / 2
    # pillonetto_bell/src/mcmc.m:73
    ## Initialize steps
    k=2 # pillonetto_bell/src/mcmc.m:76
    rejected=0 # pillonetto_bell/src/mcmc.m:77
    progress=struct() # pillonetto_bell/src/mcmc.m:78
    ## Loop
    while k < K:

        if logical_not(rem(k,500)):
            k
        # Restricted to be nonnegative
        thetas[:,k]=rmvnrnd(thetas[:,k - 1], dot ( alph,Vhat_),1,cat ([- 1,0 ],[0,- 1] ),cat(0 ,0 ).T) # pillonetto_bell/src/mcmc.m:88
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