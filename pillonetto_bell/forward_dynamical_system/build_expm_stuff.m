function [ANhat,dANhat_dqM] = build_expm_stuff(AN,tau,M)
    N = size(AN,1)-1;
    MSpl = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    
    
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

    dAN_dqM = zeros(N+1,N+1,M+2);
    % First page is q2, all zeros.
    for k=1:M+1
        dAN_dqM(:,:,k+1) = -MSpl\(N^2*(diag(-si_offdiag(k,:),-1)+diag(si_diag(k,:))+diag(-si_offdiag(k,:),1)));
    end

    
    
    ANhat_and_dANhatdqM = zeros(2*(N+1),2*(N+1),M+2);
    for k=1:M+2
        ANhat_and_dANhatdqM(:,:,k) = expm([tau*AN,tau*dAN_dqM(:,:,k);zeros(N+1),tau*AN]);
    end

    % All the ANhats should be the same -- we're calculating ANhat 2(M+2) times
    ANhat = mean(ANhat_and_dANhatdqM(1:N+1,1:N+1,:),3);
%    ANhat = mean(ANhat_and_dANhatdqM((N+2):end,(N+2):end,:),3);
    dANhat_dqM = ANhat_and_dANhatdqM(1:(N+1),(N+2):end,:);
end