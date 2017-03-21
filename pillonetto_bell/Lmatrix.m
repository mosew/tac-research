function Lmatrix = Lmatrix(theta,P,T,rkhs_eigenfile,n,tau)
    % OUTPUT:
    % P x n
    % 
    % Lmatrix(j,i) = L_i(theta,phi_j), starting everything at t=0.

    [~,eifs] = get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile);
    Lmatrix = convolve_eifs(eifs,theta,P,n,tau);
    
end