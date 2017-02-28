function p_ygivenq = p_y_given_q(y,q)
    p_ygivenq = det(2*pi*Vy_given_q(y,q))^-.5 * exp(-0.5 * y' * (v_covariance(y,q)\y));
end