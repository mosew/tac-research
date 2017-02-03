total=1:11;
p=3;
[t,u_total,y_total] = prepare_data(t_TAC_5122(total),t_BrAC_5122(total),data_corrected_TAC_5122(total),data_BrAC_5122(total),5,tau);
plot(u_total(p,:))
hold on
u_total = max(sgolayfilt(u_total',2,11),0)';
plot(u_total(p,:))