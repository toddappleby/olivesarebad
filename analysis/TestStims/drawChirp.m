function result = drawChirp(drawStruct)


bgIntensity = .5;
stepContrast = drawStruct.stepContrast;
contrastMax= drawStruct.contrastMax;
contrastMin= drawStruct.contrastMin;
preTime = drawStruct.preTime;
stepTime= drawStruct.stepTime;
tailTime = drawStruct.tailTime;
interTime = drawStruct.interTime;
frequencyTime = drawStruct.frequencyTime;
contrastTime = drawStruct.contrastTime;
contrastFreq = drawStruct.contrastFrequency;
sampleRate = drawStruct.sampleRate;
% these are all in milliseconds
timePoints = @(tMS)tMS / 1e3 * sampleRate;

totalTime = preTime + stepTime + interTime + stepTime + interTime + frequencyTime + interTime + contrastTime + interTime + tailTime;

frequencyMin = drawStruct.frequencyMin;
frequencyMax = drawStruct.frequencyMax;


%pretime
result(1,1:timePoints(preTime))=0;
%steptime
stepStimPos(1,1:timePoints(stepTime)) = stepContrast;
interStim(1,1:timePoints(interTime)) = 0;
stepStimNeg(1,1:timePoints(stepTime)) = -stepContrast;
result = [result stepStimPos interStim stepStimNeg interStim];

frequencyDelta = (frequencyMax - frequencyMin)/(frequencyTime*10)/2; %time*10 for pts--to match data?
%a frame every 16 milliseconds...but this doesn't matter when we only want
%the sweep ( and not a re-creation of the stimulus, which would require
%frame timings anyway! )  question of how this manifests in the stim tho..
%probably not diff enough to matter, tho more low peaks likely at high
%freqs of course
resultFsweep = sin(2*pi*(frequencyMin+frequencyDelta*((1:1:timePoints(frequencyTime)).^2)*1e-4));
% figure(50)
% plot(resultFsweep)

% resultFsweep = sin(2*pi*(frequencyMin+(1:1:timePoints(frequencyTime))*frequencyDelta).*(1:1:timePoints(frequencyTime))*1e-4);

%THIS WORKS, BUT IN CHIRP GEN THE DELTA IS NOT MULTIPLIED BY SECONDS?????
%also weird because these are time points into seconds but... you wouldn't
%get anywhere close to this many frames presented at 60 Hz frame rate..?
%then again, maybe you would for led..? ?????????????

%just figured it out--scaling freq Delta differently between scripts
%(smaller for the non-square version).  
% resultFsweep = sin(2*pi*((1:1:timePoints(frequencyTime))*1e-4).*(frequencyDelta*((1:1:timePoints(frequencyTime))*1e-4))); 
% resultFsweep = .5 * resultFsweep + .5;


result = [result resultFsweep];

contrastDelta = (contrastMax - contrastMin)/(contrastTime*1e-3);
resultCsweep = (contrastMin+((1:1:timePoints(contrastTime))*1e-4)*contrastDelta) .* sin(2*pi*((1:1:timePoints(contrastTime))*1e-4)*contrastFreq);

tailResult(1,1:timePoints(tailTime))=0;

result = [result interStim resultCsweep tailResult];

result = bgIntensity *result + bgIntensity;
saveResult = result';
save('chirpStim.mat','saveResult')
% figure(20)
% plot(result)
end