function EV = EV_aP_given_theta_y(theta,y,rkhs_eigenfile,P,T,n,tau,data_path)
    % INPUTS:
    % theta is a vector of parameters
    % y is observed data
    % rkhs_eigenfile is a file with output [eiv,eif] (each a single thing)
    % P is number of eigenfunctions to use
    % T is max time
    %
    % OUTPUTs:
    % EV is a cell array with 2 elements. First is a vector E, second a matrix V
    
    [eivs,~]=get_kernel_eigenstuff(P,T,rkhs_eigenfile);
    eivs = evaluate_eivs(eivs,theta);
    Vy_th_=Vy_th(theta,P,T,n,tau,rkhs_eigenfile,data_path);
    
        
    % Not sure about the time vector on this one.
    LMatrix=Lmatrix(theta,P,T,rkhs_eigenfile,n,tau)';
    %LMatrix=sample_eigenfunctions(eifs,(1/n):(1/n):T);
    
    EV = cell(1,2);
    EV{1} = diag(eivs)*(LMatrix'/Vy_th_)*y';
    EV{2} = diag(eivs)-diag(eivs)*LMatrix'*(Vy_th_\LMatrix)*diag(eivs);
    EV{2} = (EV{2} + EV{2}') / 2;
    
end