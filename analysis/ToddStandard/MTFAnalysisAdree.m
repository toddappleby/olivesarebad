
% by TA for large cell analysis
function [result,sender,spikeDetectingSTD,finalDiscardSectns] = MTFAnalysisAdree(node,params)

% epoch times(s & ms) and points
preTimeMS = node{1}.epochList.elements(1).protocolSettings.get('preTime');
preTime = preTimeMS / 1000; % now in seconds
stmTimeMS = node{1}.epochList.elements(1).protocolSettings.get('stimTime');
stmTime = stmTimeMS / 1000; %now in seconds
tailTimeMS = node{1}.epochList.elements(1).protocolSettings.get('tailTime');
tailTime = tailTimeMS / 1000; %now in seconds
samplingInterval = node{1}.epochList.elements(1).protocolSettings.get('sampleRate');
samplingInterval = 1/samplingInterval; %now in seconds/cycle
prePts = preTime / samplingInterval; %samples per length of time
stmPts = stmTime / samplingInterval;
tailPts = tailTime / samplingInterval;
allTime = preTime+stmTime+tailTime; %seconds
allPts = prePts + stmPts + tailPts;

% for dividing up the epoch for spike detection
stimTime = [preTime/samplingInterval (preTime + stmTime)/samplingInterval];

desiredSTD = params.initialSTD;

sampleRate = node{1}.epochList.elements(1).protocolSettings.get('sampleRate');


MTFMatrix = zeros(length(node),allPts);
radii = zeros(length(node));
spikeDetectingSTD = zeros(1,size(nodeData,1));
finalDiscardSectns = zeros(size(nodeData,1),3);

sender = struct();

for curNode = 1:length(node)
    
    nodeData = getSelectedData(node{curNode}.epochList,params.Amp);
    radii(curNode) = node{curNode}.epochList.elements(1).protocolSettings.get('radius');
    nodeMatrix = zeros(size(nodeData,1),size(nodeData,2));
    
    %struct params
    sender(curNode).params.sampleRate = sampleRate;
    sender(curNode).params.preTime = preTime*1e-3;
    sender(curNode).params.stimStart  = preTime*sampleRate+1;
    sender(curNode).params.stimEnd = (preTime+stmTime)*sampleRate;
    sender(curNode).params.spatialFrequency = radii(curNode);
    sender(curNode).params.temporalFrequency = node{curNode}.epochList.elements(1).protocolSettings.get('temporalFrequency');
    size(nodeData)
    for curEpoch = 1:size(nodeData,1)
        [spikes, finalSTD, finalDiscard] = convertSpikes(nodeData(curEpoch,:), ...
            stimTime, desiredSTD);
        size(spikes)
        finalDiscardSectns(curEpoch,:) = finalDiscard;
        spikeDetectingSTD(curEpoch) = finalSTD;
        nodeMatrix(curEpoch,:) = spikes;
        MTFMatrix(curNode,:) = mean(nodeMatrix,1);       
    end
    sender(curNode).data = MTFMatrix(curNode,:);   
end

    result = sMTFAnalysis(sender,1e4,'extracellular','avg');
   figure(1)
    hold on
    plot(result.uniqueSF,result.avgF1*1e4,params.color);
    title('F1 Response to Spot')
    xlabel('Spot Radius (micron)')
    ylabel('F1 Amplitude (Hz)')