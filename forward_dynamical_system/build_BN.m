function BN = build_BN(q2,N)
    MSpl = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    BN = MSpl\[zeros(N,1);q2];
    assert(all(size(BN)==[N+1,1]));
end