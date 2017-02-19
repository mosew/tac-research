% plot errors

% load('data_arrays_tau5_noeps39.mat')
% global tau
%     tau=5;

global u_total

figure
for i = 1:9
    subplot(2,5,i)
    hold on
    %for p = 1:5
        % paradigm 5
        a=cell2mat(full_deconvolved_BrACs(1:end,1,i))';
        plot(a);
    %end
    
    plot(u_total(i,:),'.')
    xlim([0 100])
    
    if i==9
        h=legend('P=4','P=8','P=16','P=32','P=64');
        pos=get(h,'position');
        set(h,'position',[0.8,0.2,pos(3:4)])
    end
end
