% test regularization parameters

reprocess_data
load('REGTEST_data_arrays_tau5_noeps39_lambda0001.mat')
h2=figure(1);
set(h2,'DefaultTextInterpreter', 'latex') 


for i = 1:9
    subplot(2,5,i)
    a=cell2mat(full_deconvolved_BrACs(3:5,i))';
    plot(a);
    hold on
    plot(u_total(i,:),'.')
    xlim([0,100])
    
    if i==9
        h=legend('P=16','P=32','P=64','interpolated BrAC');
        pos=get(h,'position');
        set(h,'position',[0.8,0.2,pos(3:4)])
    end

end
suptitle('$\lambda=10^{-4}$, $M=1$, $\tau=5$, $N=32$, paradigm 5')
