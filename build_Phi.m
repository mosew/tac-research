function Phi = build_Phi(ANhat,BNhat,us)
% CALCULATE STATE VECTORS FOR EACH EPISODE AND TIMESTEP
% Phi(:,:,i) for i=1:m are TRAINING EPISODES
% Phi(:,:,m+1) is the TEST EPISODE
    N=size(ANhat,1);
    m=size(us,1)-1;
    n=size(us,1)-1;
    Phi = zeros(N+1,n+1,m+1);
    for i=1:m+1
        Phi(:,:,i)=gen_varphi(ANhat,BNhat,us(i,:));
    end
end