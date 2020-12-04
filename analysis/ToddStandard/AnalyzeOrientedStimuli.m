function results = AnalyzeOrientedStimuli(node, params)

figure(2); clf; subplot(1, 2, 1); hold on
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

vectorX = 0;
vectorY = 0;
for CurNode = 1:length(node)
    EpochData = getSelectedData(node{CurNode}.epochList, params.Amp);
    if (params.CellAttached)
        for epoch = 1:size(EpochData, 1)
            S = spikeDetectorOnline(wavefilter(EpochData(epoch,:),6));
            SpikeTimes=S.sp;
%             [SpikeTimes, SpikeAmplitudes, RefractoryViolations] = SpikeDetection.Detector(EpochData(epoch, :));
            EpochData(epoch, :) = 0;
            EpochData(epoch, SpikeTimes) = 1/params.SamplingInterval;
        end
%     else
%         EpochData = BaselineCorrectOvation(EpochData, 1, 1000);
    end
    if (size(EpochData, 1) == 1)
        MeanResponse = filter(gausswin(params.DecimatePts), 1, EpochData)/sum(gausswin(params.DecimatePts));
    else
        MeanResponse = filter(gausswin(params.DecimatePts), 1, mean(EpochData))/sum(gausswin(params.DecimatePts));
    end
    tme = (1:length(MeanResponse))*params.SamplingInterval;
    plot(tme, MeanResponse + params.plotOffset*(CurNode-1), 'k', 'LineWidth', 2);
    pause(1);
    orientation(CurNode) = node{CurNode}.epochList.elements(1).protocolSettings.get('orientation');
    Resp(CurNode) = mean(MeanResponse);
    if (CurNode == 1)
        AveResponse = MeanResponse / max(abs(MeanResponse));
    else
        AveResponse = AveResponse + MeanResponse / max(abs(MeanResponse));
    end
    vectorX = vectorX + Resp(CurNode) * sin(2*pi*orientation(CurNode)/360);
    vectorY = vectorY + Resp(CurNode) * cos(2*pi*orientation(CurNode)/360);
end
xlabel('sec');
ylabel('sps/s');
AveResponse = AveResponse / length(node);
vectorX = vectorX / length(node);
vectorY = vectorY / length(node);
results.DSI = sqrt(vectorX.^2 + vectorY.^2) / mean(Resp);

subplot(1, 2, 2);
[sortedVal, sortedIndices] = sort(orientation);
polar(orientation(sortedIndices) * 2 * pi / 360, Resp(sortedIndices));

if (params.SaveGraphs)
    fileName = strcat(strcat(titleString), '-');
    fileName = strcat(fileName, params.fileName);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
    saveas(gcf, strcat(fileName, '.pdf'));
end



