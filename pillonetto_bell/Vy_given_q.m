function Vy_given_q = Vy_given_q(q,T,P,varargin)
    % INPUTS:
    % q is our parameter vector, [q2 q1M], which is an M+2-vector. Also a
    %   theta for the reproducing kernel?
    % T is the total length of an episode, in minutes
    % P is number of eigenfunctions to use
    % varargin, if present, is the name of a function which outputs
    %    eigenvalues/functions of a reproducing kernel of the form
    %    [eiv,eif],
    %       where eiv is a single eigenvalue and eif is a function handle
    
    data_path = '022717_234_splhr_arrays.mat';
    
    
    global tau
    load(data_path,'CNhat');
    q2=q(1);
    qM=q(2:end);
    t = 0:tau:T;
    
    % RKHS dependence
    [eivs,eifs] = get_kernel_eigenstuff(P,T,'green2_eigen');
    sampled_eifs = sample_eigenfunctions(eifs,t,P);
    
    % Data dependence
    meas_noise = Vy_given_qu(data_path);
    
    
    convolved_sampled_eifs = zeros(size(sampled_eifs));
    
    F=forward_system(q2,qM,sampled_eifs(2,:))
    
    for i = 1:P
        convolved_sampled_eifs(i,:) = CNhat*forward_system(q2,qM,sampled_eifs(i,:));
    end
    
    n=floor(T/tau);
    Vy_given_q = zeros(n,n);
    for i=1:n
        for k=1:n
            Vy_given_q(i,k) = sum(eivs(q).*convolved_sampled_eifs(:,i).*convolved_sampled_eifs(:,k)) + meas_noise(i,k);
        end
    end
    
end
    
