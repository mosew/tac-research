function p_y_given_th = p_y_given_th(y,th)
    p_y_given_th = det(2*pi*Vy_given_th(th,P))^-.5 * exp(-0.5 * y' * (Vy_given_th(th,P)\y));
end