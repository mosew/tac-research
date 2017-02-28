% Forward System
function Phi = forward_system(q2,qM,u)
    Kq = build_Kq(qM);
    AN = build_AN(Kq,q2);
    BN = build_BN(q2);
    [ANhat,~] = build_expm_stuff(AN);
    BNhat = build_BNhat(AN,ANhat,BN);
    Phi = build_Phi(ANhat,BNhat,u);
end
