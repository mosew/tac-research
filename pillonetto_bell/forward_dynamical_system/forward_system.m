% Forward System

function Phi = forward_system(q2,qM,tau,M,N,n,u)
    Kq = build_Kq(qM,N);
    AN = build_AN(Kq,q2,N);
    BN = build_BN(q2,N);
    [ANhat,~] = build_expm_stuff(AN,tau,M);
    BNhat = build_BNhat(AN,ANhat,BN);
    Phi = build_Phi(ANhat,BNhat,N,n,u);
end