function Kq = build_Kq(q1M,N)
    % This matrix gives [int_0^1 q^M_1(x) psi'^N_i(x) psi'N_j(x) ]_{ij}
    
    M = length(q1M)-1;
    
    SplineHandles = cell(1,M+1);
    for k=0:M
        SplineHandles{k+1} = @(x) (x>(k-1)/M).*(x<=k/M).*(M*x-(k-1)) + (x>k/M).*(x<(k+1)/M).*(-M*x+k+1);
    end
    si_diag=zeros(M+1,N+1); % at (k,j) gives integral of psi_k over support of psi_j
    si_offdiag =zeros(M+1,N);
    for k=1:(M+1)
        for j=1:N
            si_diag(k,j)=integral(SplineHandles{k},max(j-2,0)/N,min(j,N)/N);
            si_offdiag(k,j)=integral(SplineHandles{k},max(j-1,0)/N,j/N);
        end
        si_diag(k,N+1)=integral(SplineHandles{k},(N-1)/N,1);
    end
    
    Kq = N^2*gallery('tridiag',-q1M*si_offdiag,q1M*si_diag,-q1M*si_offdiag);
end