function Vy_th = Vy_th(theta,tau,T,P,rkhs_eigenfile,data_path)
    % INPUTS:
    % theta is a vector of parameters
    % T is the total length of an episode (in minutes?)
    % P is number of eigenfunctions to use
    % rkhs_eigenfile is a file 1x2 cell array output
    % data_path is the path of the data, for measurement noise?
    %
    % OUTPUTS:
    % Vy_given_theta is an n x n matrix, where n = length(y) = T/tau
    
    n=T/tau;    
    
    % RKHS dependence
    [eivs,eifs] = get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile);
    convolved_eifs = convolve_eifs(eifs,theta,P,n,tau);

    % Data dependence
    load(data_path,'u_total');
    meas_noise = Vy_thu(theta,n);
    
       
    Vy_th = zeros(n);
    
    for i=1:(n)
        for k=1:(n)
            Vy_th(i,k) = sum(eivs'*convolved_eifs(:,i).*convolved_eifs(:,k)) + meas_noise(i,k);
        end
    end
    
    
end