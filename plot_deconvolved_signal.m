%plot errors

reprocess_data


for i = 1:5

    
    subplot(2,3,i)
    hold on
    %m,p,la,para,i
    a=cell2mat(full_deconvolved_BrACs(1,1,1,:,i));

%       size(a)
%     
%     a=permute(a,[2,3,1]);
     a=permute(a,[2,4,1,3]);
     
    % Plot estimated BrAC
    plot(0:5:(size(a,1)*5-5),a);%['#AD1FFF','#C157FF','#D894FF','#ECCCFF']);
    
    % Plot sampled BrAC
    plot(t_BrAC_5122{i}(2:end),data_BrAC_5122{i}(2:end),'kx');
    
    % Plot interpolated BrAC
    plot(5:5:(5*length(u_total)),[u_total(i,2:end),0],'--k')
    
%     % Plot smoothed/zeroed TAC
%     plot(t_TAC_5122{i},data_TAC_5122{i},':')
    
    % Plot Gary's model's deconvolution
    %                  (training, test)
    % Note this 7 is our 6, since we omitted episode 3
    rosenEstBrACTAC = rosencode(7,i);
    pad = 6;
    plot(pad*5+rosenEstBrACTAC(:,1)*60,rosenEstBrACTAC(:,2),'-.')
    
    xlim([0 500])
    ylim([0,75])
    
    if i==5
        h=legend('trained on ep 6','trained on eps 6:7','trained on eps 6:8','trained on eps 6:9','sampled BrAC','interpolated BrAC','Est. BrAC, prev. method, trained on episode 6');
        pos=get(h,'position');
        set(h,'position',[0.55,0.2,pos(3:4)])
    end
    
end
set(0,'defaulttextinterpreter','latex')
suptitle(['$\lambda=0.1$, constant diffusivity, 2 splines/hr', char(10), '\qquad regularization on $|u|$ and $|u"|$'])