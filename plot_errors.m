% plot errors
figure
for i = 1:3
    subplot(1,3,i)
    hold on
    %for p = 1:5
        % Take M=1, paradigm 3
        a=cell2mat(full_deconvolved_BrACs(2,3,3,i));
        plot(a);
    %end
    
    plot(u_total(i,:),'.')
    
    if i==11
        h=legend('P=4','P=8','P=16','P=32','P=64');
        pos=get(h,'position');
        set(h,'position',[0.8,0.2,pos(3:4)])
    end
end
