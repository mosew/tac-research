function [Completed] = Plot_Episodes(BrAC,TAC,Episodes,E_Episodes,Participant)

Number_of_Episodes = length(Episodes);
[y_val_I,I_I] = max(Episodes(1).TAC);
Time_Displacement_I = [0];

figure(1)

plot(Episodes(1).Time,Episodes(1).TAC,'k--',Episodes(1).Time,Episodes(1).TAC,'b+')
hold on

for i = 2:Number_of_Episodes
    
    Time_Displacement_i = (Episodes(i).Start_Time - Episodes(1).Start_Time)*24;

    plot(Episodes(i).Time + Time_Displacement_i,Episodes(i).TAC,'k--',...
        Episodes(i).Time + Time_Displacement_i,Episodes(i).TAC,'b+')
     hold on
    
    [y_val_i,I_i] = max(Episodes(i).TAC);
    y_val_I = [y_val_I,y_val_i];
    I_I = [I_I,I_i];
    
    Time_Displacement_I = [Time_Displacement_I,Time_Displacement_i];

end

y_val_max = max(y_val_I);

for i = 1:Number_of_Episodes
    
    text(Episodes(i).Time(I_I(i))+ Time_Displacement_I(i),...
        1.3*y_val_max,num2str(i))
    
end

title(['Identified Drinking Episodes: Participant: ',Participant])
xlabel('Elapsed Time (hours)')
ylabel('Transdermal Alcohol level (% alcohol)')
legend (['Episodes 1 - ',num2str(Number_of_Episodes)])
axis([0,1.05*(Episodes(Number_of_Episodes).Time(end)+ Time_Displacement_I(Number_of_Episodes))...
    ,0,1.5*y_val_max])

hold off

for i = 1:Number_of_Episodes
    
    z_start = max(find((E_Episodes(i).M_TAC > .005),1),2) - 1;
    
    if isempty(z_start) 
        z_start = 1;
    end
        
    figure(i+1)
    plot(Episodes(i).Time,Episodes(i).TAC,'+')
    hold on
    plot(E_Episodes(i).Time - E_Episodes(i).Time(z_start),E_Episodes(i).M_TAC,'k-.')
    plot(E_Episodes(i).Time - E_Episodes(i).Time(z_start),E_Episodes(i).E_BrAC,'r--')
    
    xlabel('Time (hours)')
    ylabel('Alcohol Level (% alcohol)')
    title(['Participant ',Participant,' Episode ',num2str(i)])
    legend('TAC Data','Modeled TAC','Estimated BrAC')
    
    YL = ylim;             
    ylim([0,YL(2)])
    hold off
     
    warning off
    
    eval(['xlswrite(''E_BrAC_',Participant,'.xlsx'',[cellstr(datestr(x2mdate(Episodes(i).Start_Time + (E_Episodes(i).Time/24))))],i);']);
    eval(['xlswrite(''E_BrAC_',Participant,'.xlsx'',[E_Episodes(i).Time - E_Episodes(i).Time(1),max(E_Episodes(i).E_BrAC,0)],i,''B1'');']);
    
    warning on
    
end

[r1_r2_h] = BrAC_Estimator_Filter_Design(BrAC,TAC);
[Est_BrAC_TAC] = BrAC_Est_0_G_1_FD(TAC,r1_r2_h);

z_start = max(find(Est_BrAC_TAC(:,2)>.01,1) - 1,1);

if isempty(z_start) 
        z_start = 1;
    end

figure(i+2)
plot(TAC(:,1),TAC(:,2),'b+',BrAC(:,1),BrAC(:,2),'ro',...
    Est_BrAC_TAC(:,1) - Est_BrAC_TAC(z_start,1),Est_BrAC_TAC(:,2),'r--',...
    Est_BrAC_TAC(:,1) - Est_BrAC_TAC(z_start,1),Est_BrAC_TAC(:,3),'k-.')

xlabel('Time (hours')
ylabel('Alcohol Level (% alcohol)')
title(['Participant ',Participant,' Training Episode'])
legend('TAC Data','BrAC Data','Estimated BrAC','Modeled TAC')

YL = ylim;             
ylim([0,YL(2)])
hold off

Completed = 1;

end

