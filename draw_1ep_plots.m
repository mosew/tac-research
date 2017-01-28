global training test total
for i=1:10
    training = i:i+1;
    test = i;
    total=i:i+1; %this needs to be a RANGE that includes the training and test indices. used for setting max time window.
    subplot(2,5,i)
    minimize
end