load('p5122.mat')
figure
% Everything here's in 5-minute timesteps, not minutes.
for i = 1:11
    subplot(2,6,i)
    plot(TAC(start_times(i):(start_times(i)+durations(i))))
end