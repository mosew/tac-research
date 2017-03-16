function Vhat = Vhat(theta,tau,P,T,n,rkhs_eigenfile,data_path)
    % This should compute a cell array, where each cell contains
    % the Cramer-Rao lower bound for our estimates of parameter i
    %
    % nTheta is the number of parameters
    % y is the output signal
    
    nTheta = length(theta);
    
    Lmatrix_ = Lmatrix(theta,P,T,rkhs_eigenfile,n,tau);
    Vy_th_=Vy_th(theta,tau,T,P,rkhs_eigenfile,data_path);
    Vy_thu_=Vy_thu(theta,n);
    dVy_th_ = dVy_th(theta,P,T,n,rkhs_eigs,Lmatrix_);
    d2Vy_th_ = d2Vy_th(theta,n,Vy_th_,Vy_thu_);
    
    info = zeros(nTheta);
    
    
    Vhat = inv(info);
    
end