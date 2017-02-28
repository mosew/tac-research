%plot errors

reprocess_data
%load('bigtest_regH012_retry_arrays.mat');

figure
for i = 1:9
    subplot(2,5,i)
    hold on
    %m,p,la,para,i
    a=cell2mat(full_deconvolved_BrACs(1,:,1,1,i)')';
        

     size(a)
%     
%     a=permute(a,[2,3,1]);
%     a=permute(a,[2,4,1,3]);
    plot(0:5:(size(a,1)*5-5),a);
    
    plot(t_BrAC_5122{i}(2:end)-5,data_BrAC_5122{i}(2:end),'kx');
    
    plot(1:5:(5*length(u_total)),[u_total(i,2:end),0],'--k')
    
    plot(t_TAC_5122{i},data_TAC_5122{i},':')

    xlim([0 500])
    ylim([0,75])
    
    if i==9
        h=legend('2spl/hr','4spl/hr','sampled BrAC','interpolated BrAC','shifted & smoothed TAC');
        pos=get(h,'position');
        set(h,'position',[0.65,0.2,pos(3:4)])
    end
    
end
set(0,'defaulttextinterpreter','latex')
suptitle(['$\lambda=0.1$, each tested on 8 other episodes', char(10), '3 splines/hr, regularization on $|u|$ and $|u"|$'])