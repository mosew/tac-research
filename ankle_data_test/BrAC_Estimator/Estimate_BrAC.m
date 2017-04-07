function [E_Episodes] = Estimate_BrAC(BrAC,TAC,Episodes)

E_Episodes = [];

Episode = struct('Time',[],'E_BrAC',[],'M_TAC',[]); 

[r1_r2_h] = BrAC_Estimator_Filter_Design(BrAC,TAC);

for i = 1:length(Episodes)
    
    [Est_BrAC_TAC] = BrAC_Est_0_G_1_FD([Episodes(i).Time,Episodes(i).TAC],r1_r2_h);
    
    Episode.Time = Est_BrAC_TAC(:,1);
    Episode.E_BrAC = Est_BrAC_TAC(:,2);
    Episode.M_TAC = Est_BrAC_TAC(:,3);
    
    E_Episodes = [E_Episodes,Episode];
    
end

end

