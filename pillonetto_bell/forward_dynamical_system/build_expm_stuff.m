function [ANhat,dANhat_dq1] = build_expm_stuff(AN,tau)
    N = size(AN,1)-1;
    MSpl = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    M=0;
    
    SplineHandles={@(s)1};
    %     SplineHandles = cell(1,M+1);
%     for k=0:M
%         SplineHandles{k+1} = @(x) (x>(k-1)/M).*(x<=k/M).*(M*x-(k-1)) + (x>k/M).*(x<(k+1)/M).*(-M*x+k+1);
%     end
    
    si_diag=zeros(M+1,N+1); % at (k,j) gives integral of psi_k over support of psi_j
    si_offdiag =zeros(M+1,N);
    for k=1:(M+1)
        for j=1:N
            si_diag(k,j)=integral(SplineHandles{k},max(j-2,0)/N,min(j,N)/N,'ArrayValued',true);
            si_offdiag(k,j)=integral(SplineHandles{k},max(j-1,0)/N,j/N,'ArrayValued',true);
        end
        si_diag(k,N+1)=integral(SplineHandles{k},(N-1)/N,1,'ArrayValued',true);
    end

    dAN_dq1 = -MSpl\(N^2*(diag(-si_offdiag(k,:),-1)+diag(si_diag(k,:))+diag(-si_offdiag(k,:),1)));

    
%     ANhat_and_dANhatdq1 = zeros(2*(N+1),2*(N+1),M+2);
    ANhat_and_dANhatdq1 = expm([tau*AN,tau*dAN_dq1;zeros(N+1),tau*AN]);

    % All the ANhats should be the same -- we're calculating ANhat 2(M+2) times
    ANhat = ANhat_and_dANhatdq1(1:N+1,1:N+1,:);
    dANhat_dq1 = ANhat_and_dANhatdq1(1:(N+1),(N+2):end);
end