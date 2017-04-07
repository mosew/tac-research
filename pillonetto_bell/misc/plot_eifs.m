for i = 1:P
    f = eifs{i};
    plot(f(t))
    hold on
end