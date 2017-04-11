function dBNhat_dqM = build_dBNhat_dqM(AN,ANhat,dANhat_dqM,BN)
    global M N dAN_dqM MSpl
    dBNhat_dqM = zeros(N+1,1,M+2); % pages index the components of qM
    
    dBNhat_dqM(:,1,1) =(ANhat-eye(N+1))*((MSpl*AN)\[zeros(N,1);1]);
    
    for k=2:M+2
        %dBNhat_dqM(:,1,k) = -dANhat_dqM(:,:,k)*BN;
        dBNhat_dqM(:,1,k) = (dANhat_dqM(:,:,k) - (ANhat - eye(N+1))*(AN\dAN_dqM(:,:,k))) *(AN\BN);
    end
end