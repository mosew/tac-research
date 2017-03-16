function [eivs,eifs]=get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile)
    % INPUTS:
    % theta is a vector of parameters (which is only relevant here if the eigenvalues depend on theta)
    % P is the number of eigenfunctions to use
    % T is max time
    % rkhs_eigenfile is the name of a file which outputs single [number,function] pair [eiv,eif] 
    %
    % OUTPUTS:
    % eivs is a 1xP array of the eigenvalues of the RK
    % eifs is a 1xP cell array of handles of eigenfunctions of the RK
    
    eivs = zeros(P,1); %numbers
    eifs = cell(P,1); %functions of t
    
    for j=1:P
        [eivs(j),eifs{j}] = feval(rkhs_eigenfile,j,theta,T);
    end
    
end