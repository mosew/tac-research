%first-differenced

JN_and_dJN_globals
for i = 1:11
    subplot(2,6,i);
    hold on
    plot(diff(y_total(i,:)))
    plot(diff(u_total(i,:)))
    xlim([1,min(100,durations(i))])
end