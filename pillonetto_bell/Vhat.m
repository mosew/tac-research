function Vhat = Vhat(y,thetaHat,tau,P,T,n,eivs,rkhs_eigenfile,data_path)
    % This should compute an ordinary nTheta x nTheta array
    % represending the C-R lower bound for 
    %
    % y_total is the output signal
    
    
    d2logpy_th_ = d2logpy_th(y,thetaHat,tau,P,T,n,eivs,rkhs_eigenfile,data_path);
    
    d2logpth_ = -infoTheta();
    
    H = d2logpy_th_ + d2logpth_;
    
    Vhat = -inv(H);
    
end