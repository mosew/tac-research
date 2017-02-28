function Vy_given_qu = Vy_given_qu(data_path)
    
    load(data_path,'actual_errors')
        
    % HERE is where to tweak parameters such as number of training episodes
    a = cell2mat(actual_errors(1,1,1,1,:));
    
    a = permute(a,[2,5,1,3,4]);
    %a = a(1:140,:);
    Vy_given_qu = cov(a');
end