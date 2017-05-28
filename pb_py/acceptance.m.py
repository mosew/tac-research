
def acceptance(thetatry=None,thetaprev=None,y=None,tau=None,T=None,P=None,n=None,rkhs_eigenfile=None,data_path=None,*args,**kwargs):
    varargin = acceptance.varargin
    nargin = acceptance.nargin

    # Computes acceptance ratio of draw
    # INPUT:
    # thetatry and thetaprev are nTheta x 1 vectors
    # 
    # OUTPUT:
    # single number
    
    assert_(length(thetatry) == length(thetaprev))
    p_thetatry=p_theta(thetatry)
# pillonetto_bell/src/acceptance.m:11
    p_thetaprev=p_theta(thetaprev)
# pillonetto_bell/src/acceptance.m:12
    num=dot(p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path),p_thetatry)
# pillonetto_bell/src/acceptance.m:13
    den=dot(p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path),p_thetaprev)
# pillonetto_bell/src/acceptance.m:14
    A=min(1,num / den)
# pillonetto_bell/src/acceptance.m:15
    return A
    
if __name__ == '__main__':
    pass
