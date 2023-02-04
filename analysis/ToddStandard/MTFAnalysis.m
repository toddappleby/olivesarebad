% by TA for large cell analysis
function [result,sender,spikeDetectingSTD,finalDiscardSectns] = MTFAnalysis(node,params)

%directory/storage stuff
temp = node{1}.parent.parent.parent.parent.splitValue;
index = find(temp == ' ');
titleString = temp(1:index);
titleString = strcat(strcat(titleString, 'TA-'), node{1}.parent.parent.parent.parent.splitValue);
title(titleString);
storageDir = strcat(params.rootDir, titleString);

if (~exist(storageDir, 'dir'))
     mkdir(storageDir);
end
cd(storageDir);


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

sampleRate = node{1}.epochList.elements(1).protocolSettings.get('sampleRate');

% for dividing up the epoch for spike detection
stimTime = [preTime/samplingInterval (preTime + stmTime)/samplingInterval];

desiredSTD = params.initialSTD;

MTFMatrix = zeros(length(node),allPts);
radii = zeros(length(node));
% spikeDetectingSTD = zeros(1,size(nodeData,1));
% finalDiscardSectns = zeros(size(nodeData,1),3);

resetEpoch = 0;

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
   
    for curEpoch = 1:size(nodeData,1)
        
%         spikes = convertSpikesAdree(nodeData(curEpoch,:));
            [spikes, finalSTD, finalDiscard] = convertSpikesAdree(nodeData(curEpoch,:), ...
            stimTime, desiredSTD);
            nodeData(curEpoch,:) = spikes;
            finalDiscardSectns(curEpoch,:) = finalDiscard;
            spikeDetectingSTD(curEpoch) = finalSTD;
            if (resetEpoch)
                %kill epoch space
                curEpoch = holdEpochSpace;
            end
            resetEpoch = 0;
        if isempty(spikes)
            disp('deleted epoch')
            
            nodeMatrix(curEpoch,:) = [];
            resetEpoch = 1;
            holdEpochSpace = curEpoch;
        else   
            nodeMatrix(curEpoch,:) = spikes;
        end
        
  
        
        MTFMatrix(curNode,:) = mean(nodeMatrix,1);       
    end
    sender(curNode).data = MTFMatrix(curNode,:);   
end

    result = sMTFAnalysis(sender,60,'spikes','avg');
    
    figure(2)                                                                
    
    plot(result.uniqueSF,result.avgF1*1e4,params.color);
    hold on 
    plot(result.uniqueSF,result.avgF2*1e4,'r');
    if strcmp(node{1}.parent.parent.splitValue,'spot')
    title('F1 & F2 Response to Spot')
    else
        title('F1&F2 Response to Annulus')
    end
    xlabel('Spot Radius (micron)')
    ylabel('F1 Amplitude (Hz)')
    legend('F1','F2');
    
    if (params.SaveGraphs)
        figure(2)
        set(gcf,'Color','w');
        export_fig 'RF_Fig.pdf' -append
    end
end
        