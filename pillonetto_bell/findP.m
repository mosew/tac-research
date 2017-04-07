K = size(thetas,2)-burnoff;
gamma = 1.25;
ths = thetas(:,1:K);

P=150
Q = floor(gamma*P)+1;

a=amps_from_th(P,K,ths,y,rkhs_eigenfile,T,n,tau,data_path);
b=amps_from_th(Q,K,ths,y,rkhs_eigenfile,T,n,tau,data_path);

[~,eifs] = get_kernel_eigenstuff(900,T,rkhs_eigenfile);
fksCoarse = cell(P,1);
fksFine = cell(Q,1);

for i = 1:Q
    if i<=P
        fksCoarse{i} = f_from_a_eifs(a(:,i)',eifs(1:P));
    end
    fksFine{i} = f_from_a_eifs(b(:,i)',eifs(1:Q));
end


while ~isPbigenough(fksCoarse,fksFine,T)
    P=Q
    a=b;
    fksCoarse = fksFine;

    Q = floor(gamma*P)+1;
    b=amps_from_th(Q,K,ths,y,rkhs_eigenfile,T,n,tau,data_path);
    for i = 1:Q
        fksFine{i} = f_from_a_eifs(b(:,i)',eifs(1:Q));
    end
    
    
    figure
    
    coarseFcns = zeros(P,n);
    for i = 1:P
        coarseFcns(i,:) = feval(fksCoarse{i},t);
    end
    plot(mean(coarseFcns,1))
    
    hold on
    fineFcns = zeros(Q,n);
    for i = 1:Q
        fineFcns(i,:) = feval(fksFine{i},t);
    end
    plot(mean(fineFcns,1))
end    
P