P = 30;
T = 1000;
tau=5;
n=T/tau;
theta = [1,1];

t=0:tau:(T-tau);
assert(length(t)==n);

rkhs_eigenfile = 'green1_eigen';
data_path = '030117_234splhr_fixedtraining_testeps15_arrays';

load(data_path,'y_total')



[eivs,eifs]=get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile);
sampled_eifs = sample_eigenfunctions(eifs,t);
ce = convolve_eifs(eifs,theta,P,n,tau);

Vy_th_ = Vy_th(theta,tau,T,P,rkhs_eigenfile,data_path);

Lmatrix_ = Lmatrix(theta,P,T,rkhs_eigenfile,n,tau);
dLmatrix_ = dLmatrix(theta,P,n,Lmatrix_);

dVy_th_ = dVy_th(theta,P,T,n,eivs,Lmatrix_);
d2Vy_th_ = d2Vy_th(theta,P,T,n,eivs,Lmatrix_,dLmatrix_);

Vh = Vhat(y_total,theta,tau,P,T,n,rkhs_eigenfile,data_path)


