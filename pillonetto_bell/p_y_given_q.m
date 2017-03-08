function p_ygivenq = p_y_given_q(y,q)
    p_ygivenq = det(2*pi*Vy_given_q(q,P))^-.5 * exp(-0.5 * y' * (Vy_given_q(q,P)\y));
end