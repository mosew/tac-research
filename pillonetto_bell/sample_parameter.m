function thk = sample_parameter(th_km1)
    alpha_vel = get_alpha_vel(); %NOT IMPLEMENTED
    Vhat=Vhat();
    thk = normrnd(th_km1,alpha_vel*Vhat);
end