% test regularization parameters
clear
reprocess_data

 global lambda lambda2

 lambdas = [0.2, 0.3, 0.5, 1];

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
    
    figure
    
    set(0,'DefaultTextInterpreter', 'latex')
    load(sprintf('REGTEST_data_arrays_tau5_noeps39_lambda_%1.4g.mat',lambdas(l)));
    for k = 1:9
        subplot(2,5,k)
        a=cell2mat(full_deconvolved_BrACs(2:5,k))';
        plot(a);
        hold on
        plot(u_total(k,:),'.')
        xlim([0,100])

%         if k==9
%             leg=legend('P=8','P=16','P=32','P=64','interpolated BrAC');
%             pos=get(leg,'position');
%             set(leg,'position',[0.73,0.2,pos(3:4)])
%         end
    end
    title = sprintf('$\\lambda=$ %1.4g, $M=1$, $\\tau=5$, $N=32$, paradigm 5',lambdas(l));
    suptitle(title)
    leg=legend('P=8','P=16','P=32','P=64','interpolated BrAC');
    pos=get(leg,'position');
    set(leg,'position',[0.73,0.2,pos(3:4)])

end