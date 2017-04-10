function BN = build_BN(q2)
    global MSpl N
    BN = MSpl\[zeros(N,1);q2];
    assert(all(size(BN)==[N+1,1]));
end