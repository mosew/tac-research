@function
def get_kernel_eigenstuff(self):

    # INPUTS:
    # theta is a vector of parameters (which is only relevant here if the eigenvalues depend on theta)
    # P is the number of eigenfunctions to use
    # T is max time
    # rkhs_eigenfile is the name of a file which outputs single [number,function] pair [eiv,eif]
    
    # OUTPUTS:
    # eivs is a 1xP array of the eigenvalues of the RK
    # eifs is a 1xP cell array of handles of eigenfunctions of the RK
    
    eivs=cell(P,1)
# pillonetto_bell/src/get_kernel_eigenstuff.m:12
    
    eifs=cell(P,1)
# pillonetto_bell/src/get_kernel_eigenstuff.m:13
    
    
    for j in arange(1,P).reshape(-1):
        eivs[j],eifs[j]=feval(rkhs_eigenfile,j,T,nargout=2)
# pillonetto_bell/src/get_kernel_eigenstuff.m:16
    
    
    return eivs,eifs
    
if __name__ == '__main__':
    pass
    