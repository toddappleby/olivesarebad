% This function prompts the user to choose whether to keep all parts of the
% epoch (time in relation to stimulus presentation)

% INPUT  -  previousDiscardVector, a vector
% OUTPUT -  responseDiscardVector, a vector

% asa: 25/02/19

function responseDiscardVector = getUserEnteredSectionDiscard(previousDiscardVector)

discardResponse = ...
    input(['\n\nTo discard the spikes of a section, type "pre", "stim",' ...
        'or\n "post" (NO QUOTES) or if you want to restore a discarded' ...
        ' one,\n type "keep pre", "keep stim", or "keep post" (NO ' ...
        'QUOTES). \n\n Type here: '], 's');

responseDiscardVector = previousDiscardVector;

if strcmp(discardResponse, 'pre') 
    responseDiscardVector(1) = true;
elseif strcmp(discardResponse, 'keep pre')
    responseDiscardVector(1) = false;
elseif strcmp(discardResponse, 'stim') 
    responseDiscardVector(2) = true;
elseif strcmp(discardResponse, 'keep stim')
    responseDiscardVector(2) = false;
elseif strcmp(discardResponse, 'post')  
    responseDiscardVector(3) = true;
elseif strcmp(discardResponse, 'keep post')
    responseDiscardVector(3) = false;
end

end