function [ANhat,dANhat_dqM] = build_expm_stuff(AN)
    global dAN_dqM tau M N
    ANhat_and_dANhatdqM = zeros(2*(N+1),2*(N+1),M+2);
    for k=1:M+2
        ANhat_and_dANhatdqM(:,:,k) = expm([tau*AN,tau*dAN_dqM(:,:,k);zeros(N+1),tau*AN]);
    end

    % All the ANhats should be the same -- we're calculating ANhat 2(M+2) times
    ANhat = mean(ANhat_and_dANhatdqM(1:N+1,1:N+1,:),3);
%    ANhat = mean(ANhat_and_dANhatdqM((N+2):end,(N+2):end,:),3);
    dANhat_dqM = ANhat_and_dANhatdqM(1:(N+1),(N+2):end,:);
end