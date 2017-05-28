# Autogenerated with SMOP 0.32-7-gcce8558
from smop.core import *
# 

    
@function
def deiv(theta=None,P=None,T=None,*args,**kwargs):
    varargin = deiv.varargin
    nargin = deiv.nargin

    # This is for Example 7.2 of Pillonetto-Bell
    
    # OUTPUT:
    # P x nTheta array
    #   each row/column (j,i) is
    #       dLambda_j(theta)/d(theta_i)
    #   where Lambda_j is the j^th eigenvalue of the reproducing kernel.
    
    nTheta=length(theta)
# pillonetto_bell/src/example72_inputs/deiv.m:10
    deiv=zeros(P,nTheta)
# pillonetto_bell/src/example72_inputs/deiv.m:11
    js=arange(1,P)
# pillonetto_bell/src/example72_inputs/deiv.m:13
    deiv[:,1]=T ** 2 / (dot((js - 1),pi) + pi / 2) ** 2
# pillonetto_bell/src/example72_inputs/deiv.m:15
    deiv[:,2]=zeros(1,P)
# pillonetto_bell/src/example72_inputs/deiv.m:16
    return deiv
    
if __name__ == '__main__':
    pass
    