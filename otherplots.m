for i = 4:8
    subplot(2,3,i-3)
    plot((1:size(y_total,2))*5,y_total(i,:))
    hold on
    plot(t_TAC_5122{i},data_TAC_5122{i})
end