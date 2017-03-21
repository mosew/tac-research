function Vhat = Vhat(y_total,thetaHat,tau,P,T,n,eivs,rkhs_eigenfile,data_path)
    % This should compute a cell array, where each cell contains
    % the Cramer-Rao lower bound for our estimates of parameter i
    %
    % y is the output signal
    
    
    
    
    d2logpy_th_ = d2logpy_th(y_total,thetaHat,tau,P,T,n,eivs,rkhs_eigenfile,data_path);    
    
    d2logpth_ = infoTheta();
    
    m = size(y_total,1);
    
    H = d2logpy_th_ + m*d2logpth_;
    info = -H;
    
    Vhat = inv(info);
    
end