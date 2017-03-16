P = 30;
T = 1000;
rkhs_eigenfile = 'green1_eigen';

theta = [1,1];

tau=5;
data_path = '030117_234splhr_fixedtraining_testeps15_arrays';

n=T/tau;
t=0:tau:(T-tau);
assert(length(t)==n);

[eivs,eifs]=get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile);
sampled_eifs = sample_eigenfunctions(eifs,t);
ce = convolve_eifs(eifs,theta,P,n,tau);

Vy_th = Vy_th(theta,tau,T,P,rkhs_eigenfile,data_path);

Lmatrix = Lmatrix(theta,P,T,rkhs_eigenfile,n,tau);
dLmatrix = dLmatrix(theta,P,n,Lmatrix);

dVy_th = dVy_th(theta,P,T,n,eivs,Lmatrix);




%Vh = Vhat(theta,tau,T,P,n,rkhs_eigenfile,data_path);

