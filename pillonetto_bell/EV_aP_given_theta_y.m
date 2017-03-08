function EV = EV_aP_given_theta_y(theta,y,rkhs_eigenfile,P,T)
    % INPUT
    % theta is a vector of parameters
    % y is observed data
    % rkhs eigenfile is a file with output [eiv,eif] (each a single thing)
    % P is number of eigenfunctions to use
    % T is max time
    %
    % OUTPUT
    % EV is a cell array with 2 elements. First is a vector E, second a
    % matrix V
    
    [eivs,eifs]=get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile);
    Vy_given_theta=Vy_given_theta(theta,T,P,rkhs_eigenfile);
    
    
    n=length(y);
    
    % Not sure about the time vector on this one.
    LMatrix=sample_eigenfunctions(eifs,0:1/(n-1):T,P);
    
    EV = cell(1,2);
    EV{1} = diag(eivs)*LMatrix'*Vy_given_theta\y;
    EV{2} = diag(eivs)-diag(eivs)*LMatrix'*Vy_given_theta\LMatrix*diag(eivs);
    
end