display(' ');
Participant_ID = input('Enter participant ID number: ');
display(' ');
Threshold = input('Enter noise threshold (a nonnegative real number): ');
display(' ');
Min_Gap = input('Enter the minimum number of zero TAC time intervals between drinking episodes (a positive integer): ');
display(' ');

Participant = sprintf('%05d',Participant_ID);

Delete_Spreadsheet = input(['Would you like to delete the current version of the spreadsheet E_BrAC_',Participant,'.xlsx ? (Y/N) '],'s');

display(' ');

if strcmpi(Delete_Spreadsheet,'Y')
    
    warning off
    
    eval(['delete E_BrAC_',Participant,'.xlsx;']);
    
    warning on
    
end

[Completed] = BrAC_Estimator(Participant,Threshold,Min_Gap);

display('BrAC successfully estimated');
display(['Excel spreadsheet E_BrAC_',Participant,'.xlsx successfully created.']);
display(['MATLAB workspace Data_',Participant,'.mat successfully created.']);
display(' ');
