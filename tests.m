%% Run model for various values of M, P, testing paradigms, episodes
global training test tau u_total
global M P
global lambda

%% Call reprocess_data to define u_total, maybe other stuff?
%reprocess_data

%% Define cell array to hold test data

% (M, )P, (test paradigm, )test episode
b = cell(4,9);

%% Run tests

% % M=1 did best? i.e. diffusion coefficient linearly dependent on depth.

% For episodes 1 through 9 (excluded original episodes 3 and 9)
for i = 1:9
    
    
    % For testing paradigms 1 through 5
    % paradigm 1: training and testing on same episode
    % paradigm 2: testing on i, training on i : i+3 (wraparound)
    % paradigm 3: testing on i, training on i+1 : i+4 (wraparound)
    % paradigm 4: training on all episodes
    % paradigm 5: training on all except test episode
    for para = 5
        
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
        for Mind = 1
            
            M = 2^Mind;
            
            if M==1
                M=0;
            else
                M=M/2;
            end
            
            
            % For various values of P
            for Pind = (3:6)
                P = 2^Pind;
                
                
                % Run the minimization script and time it
                tic
                minimize_uspline
                [M,P,para,i]
                toc
                
                % Collect data from the run
                test_u = [u_total(i,:),0];
                [peak_est, peaktime_est] = max(u_star);
                [peak_act, peaktime_act] = max(test_u);
                
                % Define struct of collected data
                b{Pind-2,i} = struct('trained_parameters',{[q2_star,q1M_star]},...
                                             'full_deconvolved_BrAC',{u_star},...
                                             'actual_error',{u_star-test_u},...
                                             'L2_error',{sum((u_star-test_u).^2)},...
                                             'Linf_error',{max(abs(u_star-test_u))},...
                                             'AUC_error',{sum(u_star)-sum(test_u)},...
                                             'peak_time_error',{tau*(peaktime_est-peaktime_act)},...
                                             'peak_height_error',{peak_est-peak_act},...
                                             'badscale',{badscale});
                                

            end    
        end
    end
end
['Saving test results cell array']
name=sprintf('REGTEST_tau5_noeps39_lambda_%1.4g.mat',lambda);
save(name,'b')
