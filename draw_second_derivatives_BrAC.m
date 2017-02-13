% Draw second derivatives
global n
d2 = zeros(11,n-2);

for i=1:11
    subplot(2,6,i)
    plot(diff(diff(u_total(i,:))));
    d2(i,:) = diff(diff(u_total(i,:)));
    d1(i,:) = diff(u_total(i,:));

end

d2 = [d2, zeros(11,2)];
d1 = [d1, zeros(11,1)];
figure
fftd2=fft(d1');
realcoefs = abs(fftd2(1:(n/2+1),:));
plot(realcoefs)