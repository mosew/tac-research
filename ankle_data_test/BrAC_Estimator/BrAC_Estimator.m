function [Completed] = BrAC_Estimator(Participant,Threshold,Min_Gap)

close all

eval(['[BrAC_',Participant,',TAC_',Participant,',Error_Flag] = Get_Training_Data(Participant);']);

if Error_Flag 
    
    Return;
    
end

eval(['[Episodes_',Participant,',Time_TAC_',Participant,',Episode_Start_',Participant,...
    ',Episode_End_',Participant,'] = Read_and_Clean_Scram(Participant,Threshold,Min_Gap);']);
eval(['[E_Episodes_',Participant,'] = Estimate_BrAC(BrAC_',Participant,',TAC_',Participant,...
    ',Episodes_',Participant,');']);
eval(['[Completed] = Plot_Episodes(BrAC_',Participant,',TAC_',Participant,',Episodes_',Participant,...
    ',E_Episodes_',Participant,',Participant);']);

eval(['save Data_',Participant,' BrAC_',Participant,' TAC_',Participant,' Episodes_',Participant,...
      ' Time_TAC_',Participant,' E_Episodes_',Participant]);
  
end

