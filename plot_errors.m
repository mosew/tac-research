% plot errors

% load('data_arrays_tau5_noeps39.mat')
% global tau
%     tau=5;

global u_total

figure
for i = 1:5
    subplot(2,3,i)
    hold on
    %for p = 1:5
        % m,p,la,para,i
        a=cell2mat(full_deconvolved_BrACs(:,4,2,3,i))';
        %a=permute(a,[2,3,1])
        plot(a);
    %end
    
    plot(u_total(i,:),'.')
    xlim([0 100])
    ylim([0,70])
    
%     if i==5
%         h=legend('P=15','P=30','P=60');
%         pos=get(h,'position');
%         set(h,'position',[0.8,0.2,pos(3:4)])
%     end
end
