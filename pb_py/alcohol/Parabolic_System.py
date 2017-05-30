class Parabolic_System(object):

    def __init__(self):

        self.N = 32
        self.tau = 5.
        self.T = 276
        self.q2 = 1.
        self.q1 = .004
        self.n = int(self.T/self.tau)

        import numpy as np
        self.u_total = np.ones((3,self.n))
        # Here the test episode should have row index 0, not m+1
        self.y_total = None
        self.BN = None


    def build_final_matrices(self):
        import numpy as np

        N = self.N
        a = np.ones(N)
        d = 4*np.ones(N-1)
        d = np.concatenate(([2],d))
        d = np.concatenate((d,[2]))

        self.M = np.array( 1/(6*float(N)) * (np.diag(a,k=-1) + np.diag(d) + np.diag(a,k=1) ), dtype=float)

        assert(isinstance(self.M,np.ndarray))
        assert(self.M.shape == (self.N+1,self.N+1))

        self.K = np.array(N * ( np.diag(-np.ones(N),k=-1) + np.diag(d/2) + np.diag(-np.ones(N),k=1)), dtype=float )
        assert(isinstance(self.K,np.ndarray))
        assert(self.K.shape == (self.N+1,self.N+1))

    def build_Kq(self):

        self.Kq = self.q1 * self.K

    def build_AN(self):
        import numpy as np

        L=np.diag( np.concatenate( ([1],np.zeros(self.N)) ) )
        self.AN=np.linalg.solve(- self.M,(L + self.Kq))

        assert(self.AN.shape==(self.N+1,self.N+1))

    def build_BN(self):
        import numpy as np

        self.BN=np.array( np.linalg.solve(self.M,np.concatenate( ( np.zeros(self.N),[self.q2] ) ) ) ).reshape((self.N+1,1))

    def build_BNhat(self):
        import numpy as np

        self.BNhat=np.dot(np.linalg.solve(self.AN,(self.ANhat - np.identity(self.N + 1))),self.BN)
        assert(self.BNhat.shape == (self.N+1,1))


    def build_expm_stuff(self):

        # FIGURE THIS OUT AGAIN FROM THE MATH?
        import numpy as np
        from scipy.linalg import expm

        self.dAN_dq1=np.linalg.solve(- self.M,self.K)
        
        ANhat_and_dANhatdqM=expm(self.tau * np.concatenate( (np.concatenate( (self.AN, self.dAN_dq1),axis=1), np.concatenate( (np.zeros((self.N+1,self.N+1)),self.AN), axis=1) ), axis=0 ))

        self.ANhat=ANhat_and_dANhatdqM[(self.N + 1):, (self.N + 1):]
        self.dANhat_dq1=ANhat_and_dANhatdqM[:(self.N + 1),(self.N + 1):]
        self.dANhat_dq2=np.zeros((self.N+1,self.N+1))



    def build_dBNhat_dq(self):
        import numpy as np

        outp = np.concatenate((np.zeros(self.n),[1]))

        # FIGURE THIS OUT FROM THE MATH
        self.dBNhat_dq1 = self.BNhat
        self.dBNhat_dq2 = self.BNhat/self.q2



    def gen_varphi(self,u=None):

        import numpy as np

        varphi=np.zeros((self.N + 1,self.n + 1))
        
        for j in range(1,(self.n + 1)):
            assert(self.ANhat.shape == (self.N + 1,self.N+1))
            assert(self.BNhat.shape == (self.N+1,1))
            v = np.array(varphi[:,j-1])
            a = np.dot(self.ANhat,v)
            b = u[j-1]*self.BNhat.reshape(-1)
            varphi[:,j] = np.array( a+b )

        return varphi

        
    
    def build_Phi(self,us=None):
        import numpy as np

        if us is None:
            us = np.array(self.u_total)

        # CALCULATE STATE VECTORS FOR EACH EPISODE AND TIMESTEP
        # Phi[:,:,i] for i=1:m are TRAINING EPISODES
        # Phi[:,:,0] is the TEST EPISODE

        m=us.shape[0]-1

        Phi=np.zeros((self.N + 1,self.n + 1,m + 1))

        for i in range(m + 1):
            Phi[:, :, i] = self.gen_varphi(us[i, :])

        self.Phi = Phi

    
    def define_operators(self,us=None):

        self.build_final_matrices()
        self.build_Kq()
        self.build_AN()
        self.build_BN()
        self.build_expm_stuff()
        self.build_BNhat()
        self.build_Phi(us)


if __name__ == '__main__':
    P = Parabolic_System()
    P.define_operators()
    print P.BN
