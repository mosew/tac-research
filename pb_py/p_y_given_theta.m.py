# Autogenerated with SMOP 0.32-7-gcce8558
from smop.core import *
# 

    
@function
def p_y_given_theta(y=None,th=None,tau=None,T=None,P=None,n=None,rkhs_eigenfile=None,data_path=None,*args,**kwargs):
    varargin = p_y_given_theta.varargin
    nargin = p_y_given_theta.nargin

    v=Vy_th(th,P,T,n,tau,rkhs_eigenfile,data_path)
# pillonetto_bell/src/p_y_given_theta.m:2
    p_y_given_th=mvnpdf(y,0,(v + v.T) / 2)
# pillonetto_bell/src/p_y_given_theta.m:3
    return p_y_given_th
    
if __name__ == '__main__':
    pass
    