function BNhat = build_BNhat(AN,ANhat,BN)
    global N
    BNhat = AN\(ANhat-eye(N+1))*BN;
end
