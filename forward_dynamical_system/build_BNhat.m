function BNhat = build_BNhat(AN,ANhat,BN)
    N = size(AN,1)-1;
    BNhat = AN\(ANhat-eye(N+1))*BN;
end
