function [Episodes,Time_TAC,Episode_Start,Episode_End] = Read_and_Clean_Scram(FILE,Threshold,Min_Gap)

[NUM,TXT,RAW]=xlsread(FILE);

Index_r = find(strcmp(RAW(:,1),'TAC Level')) + 1;
Index_c = find(strcmp(RAW(Index_r - 1,:),'Time'));

Time_TAC = [cell2mat(RAW(Index_r:end,Index_c)),str2double(RAW(Index_r:end,1))];

I = ~isnan(Time_TAC(:,2));

Time_TAC = Time_TAC(I,:);

J = Time_TAC(:,2) >= Threshold;

Time_TAC(:,2) = J.*Time_TAC(:,2);

done = false;
pointer = 1;

Episodes = [];

Episode = struct('Start_Time',[],'Time',[],'TAC',[]); 

Episode_Start = [];
Episode_End = [];

Last_Episode_End = 0;

while (~done)
    
    start_pointer = max(find(Time_TAC(pointer:end,2),1) + pointer - 1,2);
    end_pointer = find((Time_TAC(start_pointer:end,2)==0),1) + start_pointer - 1;
    
    if ((start_pointer - Last_Episode_End < Min_Gap) & (~isempty(Episode_Start)))
        
       start_pointer = max(Episode_Start(end),2);
        
       Episode_Start = Episode_Start(1:(end-1)); 
       Episode_End = Episode_End(1:(end-1));
       Episodes = Episodes(1:(end-1));
       
    end
    
    if isempty(end_pointer)
        
        end_pointer = length(Time_TAC(:,1));
        done = true;
        
    end
    
    pointer = end_pointer;
    
    if isempty(find(Time_TAC(pointer+1:end,2),1)) 
        
        done = true;
          
    end
    
    Episode_Start = [Episode_Start,start_pointer-1];
    Episode_End = [Episode_End,end_pointer];
    
    Last_Episode_End = Episode_End(end);
    
    Episode.Start_Time = Time_TAC(start_pointer-1,1);
    Episode.Time = (Time_TAC(start_pointer-1:end_pointer,1) - Episode.Start_Time)*24; 
    Episode.TAC = Time_TAC(start_pointer-1:end_pointer,2); 
   
    Episodes=[Episodes,Episode];
    
end

Length_Time_TAC = size(Time_TAC,1);

Episode_Number = zeros(Length_Time_TAC,1);
Number_of_Episodes = length(Episode_Start);

for i = 1:Number_of_Episodes
    
    Episode_Number(Episode_Start(i):Episode_End(i)) = i;
    
end

Time_TAC=[Episode_Number,(1:Length_Time_TAC)',Time_TAC];
Time_TAC(Episode_Start,2)=-1*Time_TAC(Episode_Start,2);
Time_TAC(Episode_End,2)=-1*abs(Time_TAC(Episode_End,2));

Time_TAC(:,3) = (Time_TAC(:,3) - Time_TAC(1,3))*24;

display(['Total number of drinking episodes found: ',num2str(length(Episode_Start))])

plot(Time_TAC(:,3),Time_TAC(:,4),'k')
title(['Identified Drinking Episodes: Participant: ',FILE])
xlabel('Elapsed Time (hours)')
ylabel('Transdermal Alcohol level (% alcohol)')
legend (['Episodes 1 - ',num2str(Number_of_Episodes)])
hold off

end

