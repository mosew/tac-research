
%don't use 3 or 9
test=5;
% number of training episodes for me (episodes 6 through 9)
para=1;
% which training episode for Gary
train =6;

plot(eval(sprintf('t_TAC_5122_%i',test)),eval(sprintf('data_TAC_5122_%i',test)))
hold on
cd('C:\Users\mose\Dropbox\research\Filter Design MW')
q_1_in=eval(sprintf('alpha_5122_%i',train));
q_2_in=eval(sprintf('beta_5122_%i',train));
q_3_in=0;
BrAC=[eval(sprintf('t_BrAC_5122_%i',train))', eval(sprintf('data_BrAC_5122_%i',train))];
TAC=[eval(sprintf('t_TAC_5122_%i',train))', eval(sprintf('data_TAC_5122_%i',train))];

[q_1_out,q_2_out,q_3_out,r1_r2_h,Est_TAC_q_in,Est_TAC_q_out] = BrAC_Estimator_Filter_Design_2(q_1_in,q_2_in,q_3_in,BrAC,TAC);


plot(Est_TAC_q_in(:,1),Est_TAC_q_in(:,2)/1000)
plot(Est_TAC_q_out(:,1),Est_TAC_q_out(:,2)/1000)






cd('C:\Users\mose\Dropbox\research\matlab\me\pillonetto_bell\')
qs = permute(trained_parameters,[2,4,5,1,3]);
us = permute(full_deconvolved_BrACs,[2,4,5,1,3]);

if test>3
    test=test-1;
end
if test>9
    test=test-1;
end


% P, #training eps, test ep
q = cell2mat(qs(1,para,test));
u_star = cell2mat(us(1,para,test));

n=size(u_star,2);

q2_star = q(1);
q1M_star = q(2:3);

% Feed optimal parameters and input forward through system to check it
Phi = forward_system(q2_star,q1M_star,5,1,32,n,u_star);
y_out = zeros(1,n);
for j=1:n
    y_out(j) = CNhat * Phi(:,j,1);
end
% Account for my zero padding at beginning of each episode
plot((-6:(n-7))/12,y_out)
legend('actual TAC','modeled TAC from BrAC','modeled TAC from Garys deconvolved BrAC','modeled TAC from my deconvolved BrAC')
ylim([0,75])
xlim([0,20])





