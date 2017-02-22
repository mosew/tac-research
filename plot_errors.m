% plot errors

clear
load('bigtest_arrays.mat');

figure
for i = 1:5
    subplot(2,3,i)
    hold on
    %for p = 1:5
        % m,p,la,para,i
        a=cell2mat(full_deconvolved_BrACs(:,4,2,3,i))';
        %a=permute(a,[2,3,1])
        plot(a);
        plot(t_BrAC_5122{i}/5,data_BrAC_5122{i},'x');
    %end
    
    plot(u_total(i,:),'.')
    xlim([0 100])
    ylim([0,70])
    
    
    if i==5
        h=legend('M=0','M=1','M=4','sampled BrAC''interpolated BrAC');
        pos=get(h,'position');
        set(h,'position',[0.65,0.2,pos(3:4)])
    end
end
