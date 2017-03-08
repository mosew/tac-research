function AN = build_AN(Kq,q2,N)
    % LINEAR SPLINE MATRIX
    MSpl = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    L = diag([1,zeros(1,N)]);
    
    AN=-MSpl\(L+Kq);
    %AN = -MSpl\(L+R+Kq); % different boundary condition here leads to different AN.
    assert(all(size(AN)==[N+1,N+1]));
end