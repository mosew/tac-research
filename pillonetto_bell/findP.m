K = size(thetas,2)-burnin;
gamma = 1.25;
ths = thetas(:,1:K);

P=6
Q = floor(gamma*P)+1;

a=amps_from_th(P,K,ths,y,rkhs_eigenfile,T,n,tau,data_path);
b=amps_from_th(Q,K,ths,y,rkhs_eigenfile,T,n,tau,data_path);

[~,eifs] = get_kernel_eigenstuff(300,T,rkhs_eigenfile);
fksCoarse = cell(P,1);
fksFine = cell(Q,1);

for i = 1:K
    fksCoarse{i} = f_from_a_eifs(a(:,i)',eifs(1:P));
    fksFine{i} = f_from_a_eifs(b(:,i)',eifs(1:Q));
end


while ~isPbigenough(fksCoarse,fksFine,T)
    P=Q
    a=b;
    fksCoarse = fksFine;

    Q = floor(gamma*P)+1;
    b=amps_from_th(Q,K,ths,y,rkhs_eigenfile,T,n,tau,data_path);
    for i = 1:K
        fksFine{i} = f_from_a_eifs(b(:,i)',eifs(1:Q));
    end
    
    
    figure
    hold on
    
    coarseFcns = zeros(K,n);
    for i = 1:K
        coarseFcns(i,:) = feval(fksCoarse{i},t);
    end
    plot(mean(coarseFcns,1),'b-')
    
    lo_up_mid = confidence_limits(fksCoarse);
    loCo=lo_up_mid{1};
    hiCo=lo_up_mid{2};
    plot(loCo(t),'b--')
    plot(hiCo(t),'b--')
    
    

    fineFcns = zeros(K,n);
    for i = 1:K
        fineFcns(i,:) = feval(fksFine{i},t);
    end
    plot(mean(fineFcns,1),'r-')
    lo_up_mid = confidence_limits(fksFine);
    loFi=lo_up_mid{1};
    hiFi=lo_up_mid{2};
    plot(loFi(t),'r--')
    plot(hiFi(t),'r--')
    
end    
P