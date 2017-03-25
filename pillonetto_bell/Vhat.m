function Vhat = Vhat(y,thetaHat,tau,P,T,n,eivs,rkhs_eigenfile,data_path)
    % This should compute a cell array, where each cell contains
    % the Cramer-Rao lower bound for our estimates of parameter i
    %
    % y_total is the output signal
    
    
    d2logpy_th_ = d2logpy_th(y,thetaHat,tau,P,T,n,eivs,rkhs_eigenfile,data_path);    
    
    d2logpth_ = -infoTheta();
    
    H = d2logpy_th_ + d2logpth_;
    
    Vhat = inv(-H);
    
end