P = 15;
T = 1;
tau=1/50;
n=50;
theta = [exp(1),10];

t=0:tau:(T-tau);
assert(length(t)==n);

rkhs_eigenfile = 'green1_eigen';
data_path = 'none';

% generate i/o
actual_u = pb_7p2_example_u() ;
sampled_u = actual_u(t);


% fprintf('Generating sample output y\n')
% y = zeros(1,n);
% y(1)=0;
% for i = 1:(n-1)
%     y(i+1) = L_i([5,9.22],actual_u,i,tau);
% end

y = [0,7.88095366723674e-05,0.00103125740020056,0.00425802631203408,0.0109487106742877,0.0216727342166110,0.0363299402796804,0.0542128037925836,0.0742872217889489,0.0952438390517745,0.115773735485898,0.134695370184833,0.151113977012157,0.164215096458382,0.173782968335687,0.179509850250412,0.181500188617153,0.180170029832357,0.176002392924988,0.169616154122135,0.161547784539635,0.152696662331240,0.143452966474724,0.134612782969953,0.126611577696363,0.119861421105558,0.114578590320305,0.110958828816849,0.108905374421002,0.108205189637280,0.108760246571796,0.109861127152109,0.111364541302021,0.112570195812181,0.112986485686659,0.112348166641714,0.110225955711135,0.106496590351219,0.101132319193791,0.0942716224384557,0.0861827407805335,0.0773306326961579,0.0680696708969857,0.0589579388905687,0.0502928147264237,0.0424323522743800,0.0355219726659830,0.0295923564096101,0.0246246942156294,0.0204888398013944];
y_total = y;

% data_path = '030117_234splhr_fixedtraining_testeps15_arrays';
% 
% load(data_path,'y_total')



[eivs,eifs]=get_kernel_eigenstuff(P,T,rkhs_eigenfile);
sampled_eifs = sample_eigenfunctions(eifs,t);
ce = convolve_eifs(eifs,theta,P,n,tau);

Vy_th_ = Vy_th(theta,P,T,n,tau,rkhs_eigenfile,data_path);

Lmatrix_ = Lmatrix(theta,P,T,rkhs_eigenfile,n,tau);
dLmatrix_ = dLmatrix(theta,P,n,Lmatrix_);
d2Lmatrix_ = d2Lmatrix(theta,P,n,Lmatrix_);

dVy_th_ = dVy_th(theta,P,T,n,eivs,Lmatrix_);
d2Vy_th_ = d2Vy_th(theta,P,T,n,eivs,Lmatrix_,dLmatrix_);

d2logpy_th_ = d2logpy_th(y,theta,tau,P,T,n,eivs,rkhs_eigenfile,data_path);

Vh = Vhat(y,theta,tau,P,T,n,eivs,rkhs_eigenfile,data_path);
% eigs(Vh)
% eigs(diag(diag(Vh)))
