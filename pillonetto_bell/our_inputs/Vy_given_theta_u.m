function Vy_given_theta_u = Vy_given_theta_u(data_path)
    
    load(data_path,'actual_errors','u_total')
    % HERE is where to tweak parameters such as number of training episodes
    a = cell2mat(actual_errors(1,1,1,1,:));
    
    a = permute(a,[2,5,1,3,4]);

%     brac=trim_BrAC(u_total);
%     [~,n]=size(brac);
%     a=a(1:n,:);
    
    Vy_given_theta_u = cov(a');
    
end