% test regularization parameters
clear
reprocess_data

global lambda lambda2

lambdas = [0 10^(-4) 10^(-3) 10^(-2) 0.1 0.15];

for l = 1:length(lambdas)

%     % Everything tested with M=1, paradigm 5 only
% 
%     lambda = lambdas(l);
%     lambda2= lambdas(l);
%     ['Testing on lambda=' num2str(lambda)]
%     tests
%     test_data_analysis
%     assert(lambda==lambdas(l));
%     
% 


%load('REGTEST_data_arrays_tau5_noeps39_lambda0001.mat')
h2=figure(1);
set(h2,'DefaultTextInterpreter', 'latex') 


for k = 1:9
    subplot(2,5,k)
    a=cell2mat(full_deconvolved_BrACs(2:5,k))';
    plot(a);
    hold on
    plot(u_total(k,:),'.')
    xlim([0,100])
    
    if k==9
        h=legend('P=8','P=16','P=32','P=64','interpolated BrAC');
        pos=get(h,'position');
        set(h,'position',[0.8,0.2,pos(3:4)])
    end

end
suptitle(sprintf('$\lambda=%d$, $M=1$, $\tau=5$, $N=32$, paradigm 5',lambda),'interpreter','latex')
end