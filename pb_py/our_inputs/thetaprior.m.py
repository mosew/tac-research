# Autogenerated with SMOP 0.32-7-gcce8558
from smop.core import *
# 

    
@function
def thetaprior(theta=None,*args,**kwargs):
    varargin = thetaprior.varargin
    nargin = thetaprior.nargin

    prior=normpdf(theta,0.0046,0.0006)
# pillonetto_bell/src/our_inputs/thetaprior.m:2
    return prior
    
if __name__ == '__main__':
    pass
    