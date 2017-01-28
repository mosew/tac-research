
avgtemps=zeros(1,11);
global n
devs = zeros(11,n);
for i = 1:11
    tempi = temp(start_times(i):start_times(i)+durations(i)-1);
    avgtemps(i) = mean(tempi);
    devs(i,:) = tempi - avgtemps(i);
end

