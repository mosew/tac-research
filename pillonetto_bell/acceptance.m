function A = acceptance(qtry,qlast,y)
    % Computes acceptance ratio of draw
    A = min(1, p_y_given_q(y,qtry)*qprior(qtry) / (p_y_given_q(y,qlast)*qprior(qlast)));
end