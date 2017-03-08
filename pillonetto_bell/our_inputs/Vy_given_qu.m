function Vy_given_qu = Vy_given_qu(data_path)
    
    load(data_path,'actual_errors','u_totals')
    % HERE is where to tweak parameters such as number of training episodes
    a = cell2mat(actual_errors(1,1,1,1,:));
    
    a = permute(a,[2,5,1,3,4]);

%     brac=trim_BrAC(u_total);
%     [~,n]=size(brac);
%     a=a(1:n,:);
    
    Vy_given_qu = cov(a');
    
end