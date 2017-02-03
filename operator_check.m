% operator check

qM = ones(1,M+1);
q2 = 1;
c = zeros(1,P+1);

global training_u
Kq = build_Kq(qM);
AN = build_AN(Kq,q2);
BN = build_BN(q2);
[ANhat,dANhat_dqM] = build_expm_stuff(AN);
BNhat = build_BNhat(AN,ANhat,BN);
dBNhat_dqM = build_dBNhat_dqM(AN,ANhat,dANhat_dqM,BN);
Phi = build_Phi(ANhat,BNhat,[training_u;sample_u_spline(c)]);