# Autogenerated with SMOP 0.32-7-gcce8558
from smop.core import *
# 

    
@function
def p_theta(theta=None,*args,**kwargs):
    varargin = p_theta.varargin
    nargin = p_theta.nargin

    # this is for example 7.2 from the Pillonetto-Bell paper
    # Returns number, joint prior probability for theta evaluated at theta
    
    p_theta=dot(dot((theta[1] > 0),(theta[2] > 0)),normpdf(theta[2],10,2)) / (1 - normcdf(- 5,0,1))
# pillonetto_bell/src/example72_inputs/p_theta.m:5
    
    return p_theta
    
if __name__ == '__main__':
    pass
    