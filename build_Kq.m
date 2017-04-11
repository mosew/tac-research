function Kq = build_Kq(q1M)
    % This matrix gives [int_0^1 q^M_1(x) psi'^N_i(x) psi'N_j(x) ]_{ij}
    global N si_diag si_offdiag
    Kq = N^2*gallery('tridiag',-q1M*si_offdiag,q1M*si_diag,-q1M*si_offdiag);
end