function Vy_th = Vy_th(theta,P,T,n,tau,rkhs_eigenfile,data_path)
    % INPUTS:
    % theta is a vector of parameters
    % T is the total length of an episode (in minutes?)
    % P is number of eigenfunctions to use
    % rkhs_eigenfile is a file 1x2 cell array output
    % data_path is the path of the data, for measurement noise?
    %
    % OUTPUTS:
    % Vy_given_theta is an n x n matrix, where n = length(y) = T/tau
        
    % RKHS dependence
    [eivs,eifs] = get_kernel_eigenstuff(P,T,rkhs_eigenfile);
    convolved_eifs = convolve_eifs(eifs,theta,P,n,tau);

    % Data dependence
    if data_path == 'none'
        % Example 7.2 from P-B paper
        u = pb_7p2_example_u();
        
        meas_noise = Vy_thu([exp(1),10],n,u);
        u = u(0:tau:T);
        
        
    else
        load(data_path,'u_total');
        meas_noise = Vy_thu(theta,n,u_total);
    end
    
       
    Vy_th = zeros(n);
    
    eivs = evaluate_eivs(eivs,theta);
    
    for i=1:(n)
        for k=1:i
            Vy_th(i,k) = sum(eivs.*convolved_eifs(:,i).*convolved_eifs(:,k));
            if i~=k
                Vy_th(k,i) = Vy_th(i,k);
            end
        end
    end
    Vy_th = Vy_th + meas_noise;
%    Vy_th = Vy_th*(1.05)^2; % 5% CV noise.
%    Vy_th = (Vy_th + Vy_th')/2;
    
end