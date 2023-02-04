function results = MeanSelectedNodes(node, params)

clear Resp AveResponse;

figure(3); clf; subplot(1, 2, 1); hold on
temp = node{1}.parent.parent.parent.splitValue;
index = find(temp == ' ');
titleString = temp(1:index);
titleString = strcat(strcat(titleString, '-'), node{1}.parent.parent.splitValue);
title(titleString);
storageDir = strcat(params.rootDir, titleString)
if (~exist(storageDir, 'dir'))
     mkdir(storageDir);
end
cd(storageDir);

SamplingInterval = node{1}.epochList.elements(1).protocolSettings.get('sampleRate');       % in Hz
SamplingInterval = 1/SamplingInterval;         % now in sec
PrePts = node{1}.epochList.elements(1).protocolSettings.get('preTime') / SamplingInterval / 1000;
StmPts = node{1}.epochList.elements(1).protocolSettings.get('stimTime') / SamplingInterval / 1000;

for CurNode = 1:length(node)
    %returns meanresponse of data and response amplitude (int response) for each splitValue
    %struct also contains splitvalues
    EpochData = getSelectedData(node{CurNode}.epochList, params.Amp);
    if (params.CellAttached)
        for epoch = 1:size(EpochData, 1)
%             [SpikeTimes, SpikeAmplitudes, RefractoryViolations] = SpikeDetection.Detector(EpochData(epoch, :));
                 S = spikeDetectorOnline(wavefilter(EpochData(epoch,:),6));
                 SpikeTimes=S.sp;            
                 EpochData(epoch, :) = 0;
                 EpochData(epoch, SpikeTimes) = 1/params.SamplingInterval;
        end
    else
        EpochData = BaselineCorrectOvation(EpochData, 1, 1000);
    end
    if (size(EpochData, 1) == 1)
        MeanResponse = filter(gausswin(params.DecimatePts), 1, EpochData)/sum(gausswin(params.DecimatePts));
    else
        MeanResponse = filter(gausswin(params.DecimatePts), 1, mean(EpochData))/sum(gausswin(params.DecimatePts));
    end
    tme = (1:length(MeanResponse))*params.SamplingInterval;
    results.respAmp(CurNode) = sum(MeanResponse(PrePts+1:PrePts+StmPts)) * params.SamplingInterval;
    
   
    plot(tme, MeanResponse + params.plotOffset*(CurNode-1), 'LineWidth', 2);
    
    %params.plotColors(CurNode-1) -- removed from plot above because
    %causing issues
    
    
    legendString{CurNode} = node{CurNode}.splitValue;
    results.meanResponse(CurNode, :) = MeanResponse;
    if (isnumeric(node{CurNode}.splitValue))
        results.splitValue(CurNode) = node{CurNode}.splitValue;
    else
        results.splitValue(CurNode) = CurNode;
    end    
end
xlabel('sec');
ylabel('sps/s');
if (~isnumeric(node{CurNode}.splitValue))
    legend(legendString);
end

subplot(1, 2, 2);
plot(results.splitValue, results.respAmp, 'o');
ylabel('integrated resp');
xlabel('split value');

if (params.SaveGraphs)
    fileName = strcat(strcat(titleString), '-');
    fileName = strcat(fileName, params.fileName);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9.125, 7.25], 'PaperUnits', 'Inches', 'PaperSize', [9.125, 7.25])
    saveas(gcf, strcat(fileName, '.pdf'));
end
