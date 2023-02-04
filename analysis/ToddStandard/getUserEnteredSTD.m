% This function prompts the user to  enter in an STD value for spike 
% counting a single experimental epoch. 

% INPUT  -  previousSTD, a double
% OUTPUT -  responseSTD, a double

% asa: 25/02/19

function responseSTD = getUserEnteredSTD(previousSTD)

stdResponse = ...
    input('Type in a new standard deviation for this epoch: ', 's');

if strcmp(stdResponse, '')
    responseSTD = previousSTD;
else
    responseSTD = regexprep(stdResponse,' +','');
    responseSTD = str2double(responseSTD); % cast what user typed as double
end

end