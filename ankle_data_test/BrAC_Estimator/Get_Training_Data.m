function [BrAC,TAC,Error_Flag] = Get_Training_Data(Participant)

BrAC = [];
TAC = [];
Error_Flag = false;

[NUM,TXT,RAW]=xlsread('TAC-BrAC Conversion Dataset');

Start_Row = find((NUM(:,1)== str2num(Participant)));

if isempty(Start_Row)

    display(['Error: Training data for participant ',Participant,' not found.']) 
    Error_Flag = true;
   
    return
    
end

End_Row = find((~isnan(NUM((Start_Row + 1):end,1))),1) + Start_Row - 1;

if isempty(End_Row) End_Row = length(Num(:,1));
    
end

BrAC = NUM(Start_Row:End_Row,[2,3]);
I = ~isnan(BrAC(:,1));
BrAC = BrAC(I,:);
BrAC(:,1) = (BrAC(:,1) - BrAC(1,1))*24;

TAC = NUM(Start_Row:End_Row,[4,5]);
I = ~isnan(TAC(:,1));
TAC = TAC(I,:);
TAC(:,1) = (TAC(:,1) - TAC(1,1))*24;


end

