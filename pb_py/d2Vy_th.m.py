# Autogenerated with SMOP 0.32-7-gcce8558
from smop.core import *
# 

    
@function
def d2Vy_th(theta=None,P=None,T=None,n=None,eivs=None,Lmatrix=None,dLmatrix=None,*args,**kwargs):
    varargin = d2Vy_th.varargin
    nargin = d2Vy_th.nargin

    # OUTPUT:
    # n x n x nTheta x nTheta
    
    eivs=evaluate_eivs(eivs,theta)
# pillonetto_bell/src/d2Vy_th.m:6
    nTheta=length(theta)
# pillonetto_bell/src/d2Vy_th.m:8
    deivs=deiv(theta,P,T)
# pillonetto_bell/src/d2Vy_th.m:10
    d2Lmatrix_=d2Lmatrix(theta,P,n,Lmatrix)
# pillonetto_bell/src/d2Vy_th.m:11
    d2eiv_=d2eiv(theta,P,T)
# pillonetto_bell/src/d2Vy_th.m:12
    d2Vy_th=zeros(n,n,nTheta,nTheta)
# pillonetto_bell/src/d2Vy_th.m:14
    for i in arange(1,n).reshape(-1):
        for k in arange(1,n).reshape(-1):
            for s in arange(1,nTheta).reshape(-1):
                for r in arange(1,nTheta).reshape(-1):
                    d2Vy_th[i,k,s,r]=sum(multiply(multiply(d2eiv_[:,s,r],Lmatrix[:,i]),Lmatrix[:,k]) + multiply(deivs[:,r],(multiply(dLmatrix[:,i,s],Lmatrix[:,k]) + multiply(Lmatrix[:,i],dLmatrix[:,k,s]))) + multiply(deivs[:,s],(multiply(dLmatrix[:,i,r],Lmatrix[:,k]) + multiply(Lmatrix[:,i],dLmatrix[:,k,r]))) + multiply(eivs,(multiply(d2Lmatrix_[:,i,s,r],Lmatrix[:,k]) + multiply(Lmatrix[:,i],d2Lmatrix_[:,k,s,r]) + multiply(dLmatrix[:,i,s],dLmatrix[:,k,r]) + multiply(dLmatrix[:,i,r],dLmatrix[:,k,s]))))
# pillonetto_bell/src/d2Vy_th.m:20
    
    
    d2Vy_th=d2Vy_th + d2Vy_thu(theta,P,n)
# pillonetto_bell/src/d2Vy_th.m:31
    return d2Vy_th
    
if __name__ == '__main__':
    pass
    