% plot errors
figure
for i = 1:11
    subplot(2,6,i)
    hold on
    for p = 1:5
        % Take M=1, paradigm 2
        a=cell2mat(actual_errors(2,p,2,i));
        plot(a);
    end
    if i==11
        h=legend('P=4','P=8','P=16','P=32','P=64');
        pos=get(h,'position');
        set(h,'position',[0.8,0.2,pos(3:4)])
    end
end
