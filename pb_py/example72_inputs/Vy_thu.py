# Autogenerated with SMOP 0.32-7-gcce8558
from smop.core import *
# 

    
@function
def Vy_thu(theta=None,n=None,u=None,*args,**kwargs):
    varargin = Vy_thu.varargin
    nargin = Vy_thu.nargin

    
    tau=1 / n
# ./pillonetto_bell/src/example72_inputs/Vy_thu.m:4
    conv_u=zeros(1,n)
# ./pillonetto_bell/src/example72_inputs/Vy_thu.m:5
    for i in arange(1,n).reshape(-1):
        Vy_thu[i,i]=(dot(0.05,L_i(theta,u,i,tau))) ** 2
# ./pillonetto_bell/src/example72_inputs/Vy_thu.m:9
    
    
    Vy_thu=dot(dot(0.05,0.05),diag(conv_u) ** 2)
# ./pillonetto_bell/src/example72_inputs/Vy_thu.m:13
    return Vy_thu
    
if __name__ == '__main__':
    pass
    