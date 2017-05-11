function dAN_dq1 = build_dAN_dqM(AN)
    N = size(AN,1)-1;
    MSpl = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    
    M=0;
    
    si_diag=zeros(M+1,N+1); % at (k,j) gives integral of psi_k over support of psi_j
    si_offdiag =zeros(M+1,N);
    for j=1:N
        si_diag(j)=min(j,N)/N-max(j-2,0)/N;
        si_offdiag(j)=j/N-max(j-1,0)/N;
    end
    si_diag(N+1)=1/N;

    dAN_dq1 = -MSpl\(N^2*(diag(-si_offdiag,-1)+diag(si_diag)+diag(-si_offdiag,1)));

end