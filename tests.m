%% Run model for various values of M, P, testing paradigms, episodes
global training test tau u_total
global M P
global lambda lambda2

%% Call reprocess_data to define u_total, maybe other stuff?
reprocess_data

%% Define cell array to hold test data

% tau=5,N=32
% M, P, lambda, test paradigm, test episode
b = cell(1,1,3,1,5);

Ms=[1];
Ps=[120]; %splines per hour = P spl/ep * 1 ep/240step * 1step/5min * 60min/1hr = 60P/(240*5) = P/20.
lambdas=[0.01,0.05,0.1];

%% Run tests


% For episodes 1 through 5 (excluded original episodes 3 and 9)
testeps=1:5;
for i = testeps
    
    
    % For testing paradigms 1 through 5
    % paradigm 1: training and testing on same episode
    % paradigm 2: testing on i, training on i : i+3 (wraparound)
    % paradigm 3: testing on i, training on i+1 : i+4 (wraparound)
    % paradigm 4: training on all episodes
    % paradigm 5: training on all except test episode
    paradigms=1;
    for paraInd = paradigms
        
        para = para - para(paraInd(1)) + 1;
        
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
        for Mind = 1:length(Ms)
            M = Ms(Mind);
            
            % For various values of P
            for Pind = 1:length(Ps)
                P = Ps(Pind);
                
                % For various regularizations
                for lambdaInd = 1:length(lambdas)
                    lambda = lambdas(lambdaInd);
                    lambda2= lambdas(lambdaInd);
                
                    
                    % Run the minimization script and time it
                    tic
                    minimize_uspline
                    [M,P,lambda,para,i]
                    toc

                    % Collect data from the run
                    test_u = [u_total(i,:),0];
                    [peak_est, peaktime_est] = max(u_star);
                    [peak_act, peaktime_act] = max(test_u);
                    
                    if Mind==1
                        M=0;
                    end
                    
                    % Define struct of collected data
                    b{Mind,Pind,lambdaInd,paraInd,i} = struct('tau_M_P_lambda_paradigm_test',{[tau,M,P,lambda,para,i]},...
                                                 'trained_parameters',{[q2_star,q1M_star]},...
                                                 'full_deconvolved_BrAC',{u_star},...
                                                 'actual_error',{u_star-test_u},...
                                                 'L2_error',{sum((u_star-test_u).^2)},...
                                                 'Linf_error',{max(abs(u_star-test_u))},...
                                                 'AUC_error',{sum(u_star)-sum(test_u)},...
                                                 'peak_time_error',{tau*(peaktime_est-peaktime_act)},...
                                                 'peak_height_error',{peak_est-peak_act},...
                                                 'badscale',{badscale});
                                             
                    if Mind==1
                        if q1M_star(1)~=q1M_star(2)
                            fprintf('CONSTANT DIFFUSIVITY RUN BROKEN')
                        end
                    end
                    
                end
            end    
        end
    end
end
%['Saving test results cell array']
%save('bigtest_regH012_retry.m','b')
