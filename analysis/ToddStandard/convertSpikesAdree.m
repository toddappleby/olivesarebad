
function [response, finalSTD, finalDiscard] = ...
    convertSpikesAdree(originalResp, stimTime, originalSTD)

% [SpikeTimes, amplitudes, ~] = SpikeDetection.Detector(resp);
%Convert cell into a double matrix.

noEpochDecision = true;
workingSTD = originalSTD;
epochDivided = {[1:stimTime(1)], [stimTime(1)+1:stimTime(2)], ...
    [stimTime(2)+1:length(originalResp)]};

discardSection = [false false false]; %corresponds to sections of epochs


while noEpochDecision
    
    resp = originalResp;
    
    S = spikeDetectorOnlineAdree(resp, workingSTD);
    SpikeTimes = S.sp; 
    amplitudes = S.spikeAmps;
    
    indexEpoch = zeros(size(resp));
    amplitudeEpoch = zeros(size(resp));
    
    for n = 1:size(resp,1)
        for m = 1:size(SpikeTimes,2)
            temp = SpikeTimes(1,m);
            sample = amplitudes(1,m);
            indexEpoch(n,temp) = 1; % not sure i need the spike check below, so took it related stuff here out
            amplitudeEpoch(n,temp) = sample;
        end
    end
    
    % This lets you remove spikes from particular time intervals in the
    % epoch.
    if any(discardSection)
        sectionsToZero = find(discardSection);
        indicesToZero = [];
        for i = 1:length(sectionsToZero)
            indicesToZero = [indicesToZero epochDivided{sectionsToZero}];
        end
        
        indexEpoch(indicesToZero) = 0;
        %amplitudeEpoch(indicesToZero) = 0;
        
    end
    
    % Create vectors where 1 denotes the time of a spike. (note: did this
    % above)
    response = indexEpoch;
    
    test = response .* resp;
    
    test(test == 0) = NaN; % remove irrelevant zeros for plot
    
    figure(4) % plot individual epochs
    plot(1:size(response,2),originalResp)
    hold on
    scatter(1:size(response,2),test) % results from spike detection
    hold off
    
    
    
    % respond to epochs here
    respondedToEpoch = false;
    while ~respondedToEpoch
        pause;
        userResponse = get(gcf, 'CurrentKey');
            
        % modify the STD for spike detection in the terminal
        if strcmp(userResponse,'space')
            workingSTD = getUserEnteredSTD(workingSTD);
            respondedToEpoch = true;
            noEpochDecision = true;
            
        % get rid of a section of spikes
        elseif strcmp(userResponse, 'x')
            discardSection = getUserEnteredSectionDiscard(discardSection);
            respondedToEpoch = true;
            noEpochDecision = true;
            
        % get rid of the epoch
        elseif strcmp(userResponse,'backspace')
            workingSTD = NaN;
            discardSection = [true true true];
            response = [];
            respondedToEpoch = true;
            noEpochDecision = false;
            
        % keep epoch and STD for spike detection as is
        elseif strcmp(userResponse,'return')
            respondedToEpoch = true;
            noEpochDecision = false;
        
        
        elseif strcmp(userResponse,'q')
            disp('wtf')
        end
            
        
    end
    
end

finalSTD = workingSTD;
finalDiscard = discardSection;

end