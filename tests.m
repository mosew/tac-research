% tests

global training test total tau
total = 1:11;
global M P



% test paradigm, M, P, test episode
% paradigm 1: training and testing on same episode
% paradigm 2: training on all episodes
% paradigm 3: training on all except test episode
test_data_tau2 = cell(3,5,6,11);



for Mind = (0:4)
    for Pind = (2:7)
        for i = 1:11
            tic
            training = [1:(i-1),(i+1):11];
            test = i;
            M = 2^Mind;
            P = 2^Pind;
            if M==1
                M=0;
            else
                M=M/2;
            end
            minimize_uspline
            u_star_orig = u_star;
            test_u_orig = test_u;
            u_star_end = u_star(1:durations(i));
            test_u_end = test_u(1:durations(i));
            [peak_est, peaktime_est] = max(u_star_end);
            [peak_act, peaktime_act] = max(test_u_end);
            test_data_tau2{3,Mind+1,Pind-1,i} = struct('trained_parameters',{[q2_star,q1M_star]},...
                                         'full_deconvolved_BrAC',{u_star_orig},...
                                         'actual_error',{u_star_end-test_u_end},...
                                         'L2_error',{sum((u_star_end-test_u_end).^2)},...
                                         'Linf_error',{max(abs(u_star_end-test_u_end))},...
                                         'AUC_error',{sum(u_star_end)-sum(test_u_end)},...
                                         'peak_time_error',{tau*(peaktime_est-peaktime_act)},...
                                         'peak_height_error',{peak_est-peak_act},...
                                         'badscale',{badscale});
            [M,P,i]
            toc
        end
    end
end
