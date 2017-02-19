%% Run model for various values of M, P, testing paradigms, episodes

global training test tau u_total
global M P
tau = 5;

%% Call reprocess_data to define u_total and stuff
reprocess_data

%% Define cell array to hold test data
% (M,) P, (test paradigm,) test episode
REGTEST_test_data_tau5_noeps39_lambda01 = cell(5,1,9);

%% Run tests

% M=1 did best, i.e. diffusion coefficient linearly dependent on depth.
M=1;

% For episodes 1 through 9 (excluded original episodes 3 and 9)
for i = 1:9
    
    % For testing paradigms 1 through 5
    % paradigm 1: training and testing on same episode
    % paradigm 2: testing on i, training on i : i+3 (wraparound)
    % paradigm 3: testing on i, training on i+1 : i+4 (wraparound)
    % paradigm 4: training on all episodes
    % paradigm 5: training on all except test episode
    for para = 5:5
        if para == 1
            training = i;
        end
        if para==2
            if i<=6
                training = i:i+3;
            else
                training = [i:min(i+3,9),1:(i+3-9)];
            end
        end
        if para==3
            if i<=5
                training = (i+1):i+4;
            else
                training = [(i+1):min(i+4,9),1:(i+4-9)];
            end
        end
        if para==4
            training = 1:9;
        end
        if para == 5
            training = [1:(i-1),(i+1):9];
        end
        
        % Test episode i
        test = i;

        % For various values of M
        %for Mind = (0:4)
%            if M==1
%                 M=0;
%             else
%                 M=M/2;
%             end
%                 
%             M = 2^Mind;

            % For various values of P
            for Pind = (2:6)
                P = 2^Pind;
                
                
                % Run the minimization script and time it
                tic
                minimize_uspline
                [P,para,i]
                toc
                
                % Collect data from the run
                u_star_orig = u_star;
                test_u_orig = u_total(i,:);
                u_star_end = u_star(1:durations(i));
                test_u_end = u_total(i,1:durations(i));
                [peak_est, peaktime_est] = max(u_star_end);
                [peak_act, peaktime_act] = max(test_u_end);
                
                % Define struct of collected data
                REGTEST_test_data_tau5_noeps39_lambda01{Pind-1,1,i} = struct('trained_parameters',{[q2_star,q1M_star]},...
                                             'full_deconvolved_BrAC',{u_star_orig},...
                                             'actual_error',{u_star_end-test_u_end},...
                                             'L2_error',{sum((u_star_end-test_u_end).^2)},...
                                             'Linf_error',{max(abs(u_star_end-test_u_end))},...
                                             'AUC_error',{sum(u_star_end)-sum(test_u_end)},...
                                             'peak_time_error',{tau*(peaktime_est-peaktime_act)},...
                                             'peak_height_error',{peak_est-peak_act},...
                                             'badscale',{badscale});
                                         

            end    
        %end
    end
end
