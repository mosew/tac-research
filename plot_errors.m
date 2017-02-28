%plot errors

reprocess_data
%load('smalltest_alleps_regH02_arrays.mat');

figure
for i = 1:9
    subplot(2,5,i)
    hold on
    %m,p,la,para,i
    a=cell2mat(actual_errors(2,1,1,2,i)')';
        
    xlim([0,700])
     size(a)
%     
%     a=permute(a,[2,3,1]);
%     a=permute(a,[2,4,1,3]);
    plot(0:5:(size(a,1)*5-5),a);
    

    
%     if i==9
%         h=legend('2spl/hr','4spl/hr','sampled BrAC','interpolated BrAC','shifted & smoothed TAC');
%         pos=get(h,'position');
%         set(h,'position',[0.65,0.2,pos(3:4)])
%     end
    
end
set(0,'defaulttextinterpreter','latex')
suptitle(['$\lambda=0.1$, each tested on 8 other episodes', char(10), '3 splines/hr, regularization on $|u|$ and $|u"|$'])