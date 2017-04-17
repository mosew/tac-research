%% Run MY model for various values of M, P, testing paradigms, episodes

%% Set up globals
global training test tau u_total
global M P
global lambda lambda2

%% Call reprocess_data to define u_total, maybe other stuff?
fprintf('Processing raw data...\n')
reprocess_data

%% Define test parameters
% tau=5,N=32

% M, P, lambda, test paradigm, test episode

Ms=[0];
Ps=[40]; %splines per hour = P spl/ep * 1 ep/240step * 1step/5min * 60min/1hr = 60P/(240*5) = P/20; 1 spline per 1200/P minutes.
lambdas=[0.1];
traineps=nchoosek(1:9,4);
testeps=1:9;

%% Define cell array to hold test data
fprintf('Creating empty cell array\n')
b = cell(length(traineps),length(testeps));
 M=0;
 P=40;
lambda=0.1;
lambda2=0.1;
fprintf('Done creating empty cell array\n')
numRuns=numel(b);
thisRun=0;
rtTot=0;

%% Run tests


% Excluded original episodes 3 and 9
% Excluded 3 because (adjusted) TAC signal never returns to 0
% Excluded 9 because Susan took off the bracelet during the episode
for i = 1:9
    test=i;
    % Current test episode
    for tr = 1:length(traineps)
        training = traineps(tr,:);
    
%     % For testing paradigms 1 through 5
%     % paradigm 1: training and testing on same episode
%     % paradigm 2: testing on i, training on i : i+3 (wraparound)
%     % paradigm 3: testing on i, training on i+1 : i+4 (wraparound)
%     % paradigm 4: training on all episodes
%     % paradigm 5: training on all except test episode
%     % paradigm 6: training on episodes i+1:i+3
%     for paraInd = 1:length(paradigms)
%         para = paradigms(paraInd);
%         if para == 1
%             training = i;
%         end
%         if para==2
%             if i<=6
%                 training = i:i+3;
%             else
%                 training = [i:min(i+3,9),1:(i+3-9)];
%             end
%         end
%         if para==3
%             if i<=5
%                 training = (i+1):i+4;
%             else
%                 training = [(i+1):min(i+4,9),1:(i+4-9)];
%             end
%         end
%         if para == 4
%             training = 1:9;
%         end
%         if para == 5
%             training = [1:(i-1),(i+1):9];
%         end
%         if para == 6
%             if i<=6
%                 training = (i+1):i+3;
%             else
%                 training = [(i+1):min(i+3,9),1:(i+3-9)];
%             end
%         end
%         if para == 7
%             training = 6;
%         end
%         if para == 8
%             training = 6:7;
%         end
%         if para == 9
%             training = 6:8;
%         end
%         if para == 10
%             training = 6:9;
%         end
        
        
%         % Test episode i
%         test = i;

%         
%         % For various values of M
%         for Mind = 1:length(Ms)
%             M = Ms(Mind);
%             
%             % For various values of P
%             for Pind = 1:length(Ps)
%                 P = Ps(Pind);
%                 
%                 % For various regularizations
%                 for lambdaInd = 1:length(lambdas)
%                     lambda = lambdas(lambdaInd);
%                     lambda2= lambdas(lambdaInd);
%                 
                    
                    % Run the minimization script and time it
                    tic
                    fprintf('Testing on M=%i, P=%i, lambda=%0.3g, training episode(s) index=%i, test episode=%i \n',M,P,lambda,tr,test)
                    minimize_uspline
                    r=toc;

                    thisRun = thisRun +1;
                    fprintf('Completed %i / %i \n', thisRun,numRuns)
                    rtTot = rtTot+r;
                    fprintf('Average time per run=%f s\n',rtTot/thisRun);

                    % Collect data from the run
                    test_u = [u_total(test,:),0];
                    [peak_est, peaktime_est] = max(u_star);
                    [peak_act, peaktime_act] = max(test_u);
                    
                    % This should be commented out if we aren't testing
                    % constant diffusivity parameters.
                    if Mind==1
                        M=0;
                    end
                    
                    % Define and collect struct of collected data
                    fprintf('Filling in cell %i,%i\n\n',tr,test);
                    b{tr,test} = struct('tau',{tau},'M',{M},'P',{P},'lambda',{lambda},...
                                                 'training_episodes',{training},'test_episode',{test},...
                                                 'trained_parameters',{[q2_star,q1M_star]},...
                                                 'full_deconvolved_BrAC',{u_star},...
                                                 'actual_error',{u_star-test_u},...
                                                 'L2_error',{sum((u_star-test_u).^2)},...
                                                 'Linf_error',{max(abs(u_star-test_u))},...
                                                 'AUC_sq_error',{(sum(u_star)-sum(test_u)).^2},...
                                                 'peak_time_sq_error',{(tau*(peaktime_est-peaktime_act)).^2},...
                                                 'peak_height_sq_error',{(peak_est-peak_act).^2},...
                                                 'badscale',{badscale});
                                             
                    if Mind==1
                        if q1M_star(1)~=q1M_star(2)
                            fprintf('CONSTANT DIFFUSIVITY RESTRICTION BROKEN\n')
                        end
                    end
                    
                    
%             end
%         end    
%         end
    end
end
fprintf('Saving test results cell array\n')
save('0416_2splhr_4training_alltesteps.mat','b')