function Kq = build_Kq(q1M,N)
    % This matrix gives [int_0^1 q^M_1(x) psi'^N_i(x) psi'N_j(x) ]_{ij}
    % HARD-CODED M=0
    M = 0;
    
    si_diag = zeros(1,N+1);
    si_offdiag = zeros(1,N);
    for k=1:(M+1)
        for j=1:N
            si_diag(k,j)=min(j,N)/N-max(j-2,0)/N;
            si_offdiag(k,j)=j/N-max(j-1,0)/N;
        end
        si_diag(k,N+1)=1/N;
    end
    
    Kq = N^2*gallery('tridiag',-q1M*si_offdiag,q1M*si_diag,-q1M*si_offdiag);
end