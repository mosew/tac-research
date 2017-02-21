% test regularization parameters
clear
reprocess_data


global lambda lambda2
lambdas = [0,0.001,0.01,0.1,0.2,0.4,0.8];


% Everything tested with M=1, paradigm 5 only

for l = 1:length(lambdas)

    lambda = lambdas(l);
    lambda2= lambdas(l);
    ['Testing on lambda=' num2str(lambda)]
    tests
    test_data_analysis
    assert(lambda==lambdas(l));
    
    
    set(0,'DefaultTextInterpreter', 'latex')
    figure
    try
        load(sprintf('REGTEST_data_arrays_tau5_noeps39_lambda_%1.4g.mat',lambdas(l)));
    catch
        try
            load(sprintf('REGTEST_data_arrays_tau5_noeps39_lambda_10^%1.4g.mat',log10(lambdas(l))));
        catch
            load(sprintf('REGTEST_data_arrays_tau5_noeps39_lambda_10^%1.4g.mat',round(log10(lambdas(l)))));
        end

    end
    for k = 1:9
        subplot(2,5,k)
        a=cell2mat(full_deconvolved_BrACs(:,k))';
        plot(a);
        hold on
        plot(u_total(k,:),'.')
        xlim([0,100])
        ylim([0,70])

    end
    title = sprintf('$\\lambda=$ %1.4g, $M=1$, $\\tau=5$, $N=32$, paradigm 5',lambdas(l));
    suptitle(title)
    leg=legend('P=16','P=32','P=64','interpolated BrAC');
    pos=get(leg,'position');
    set(leg,'position',[0.73,0.2,pos(3:4)])

end