function results = AnalyzeChirpStimulus(node, params)

scrsz = get(0, 'ScreenSize');

figure(2); clf; hold on
figure(3); clf;
temp = node.parent.parent.splitValue;
index = find(temp == ' ');
titleString = temp(1:index);
titleString = strcat(strcat(titleString, '-'), node.parent.splitValue);
storageDir = strcat(params.rootDir, titleString)
if (~exist(storageDir, 'dir'))
     mkdir(storageDir);
end
cd(storageDir);

SamplingRate = 1e3/node.epochList.elements(1).protocolSettings.get('sampleRate');
FreqPts = node.epochList.elements(1).protocolSettings.get('frequencyTime') / SamplingRate;
ContrastPts = node.epochList.elements(1).protocolSettings.get('contrastTime') / SamplingRate;
InterPts = node.epochList.elements(1).protocolSettings.get('interTime') / SamplingRate;
if (params.CellAttached == 1)
    PrePts = node.epochList.elements(1).protocolSettings.get('preTime') / SamplingRate;
else
    PrePts = InterPts;
end
StepPts = node.epochList.elements(1).protocolSettings.get('stepTime') / SamplingRate;
ContrastPeriod = 1e3 / (SamplingRate * node.epochList.elements(1).protocolSettings.get('contrastFrequency'));
ContrastCycles = round(ContrastPts / ContrastPeriod);
ContrastMin = node.epochList.elements(1).protocolSettings.get('contrastMin');
ContrastMax = node.epochList.elements(1).protocolSettings.get('contrastMax');
FreqMin = node.epochList.elements(1).protocolSettings.get('frequencyMin');
FreqMax = node.epochList.elements(1).protocolSettings.get('frequencyMax');

% mean chirp response
TempEpochData = getSelectedData(node.epochList, params.Amp);

if (params.CellAttached == 0)
    EpochData = ApplyFrequencyCutoffOvation(TempEpochData, params.FrequencyCutoff, params.SamplingInterval);
else
    EpochData = TempEpochData;
end
% EpochData = BaselineCorrectOvation(EpochData, 1, PrePts);
%test
if (params.CellAttached)
    for epoch = 1:size(EpochData, 1)
        S = spikeDetectorOnline(wavefilter(EpochData(epoch,:),6));
        SpikeTimes=S.sp;
%         [SpikeTimes, SpikeAmplitudes, RefractoryViolations] = SpikeDetection.Detector(EpochData(epoch, :));
        EpochData(epoch, :) = 0;
        EpochData(epoch, SpikeTimes) = 1/params.SamplingInterval;
    end
    if (size(EpochData, 1) == 1)
        MeanResponse = filter(gausswin(params.DecimatePts), 1, EpochData)/sum(gausswin(params.DecimatePts));
    else
        MeanResponse = filter(gausswin(params.DecimatePts), 1, mean(EpochData))/sum(gausswin(params.DecimatePts));
    end
else
    if (size(EpochData, 1) == 1)
        MeanResponse = EpochData;
    else
        MeanResponse = mean(EpochData);
    end
end

figure(3);
subplot(4, 1, 1);
tme = (1:length(MeanResponse))*params.SamplingInterval;
plot(tme, MeanResponse, 'LineWidth', 2);
axis tight
pause(1);
title(titleString);
if (params.SaveToIgor)
    CurFigPanel = gca;
    makeAxisStruct(CurFigPanel,'ChirpMeanResponse','test')
end

% components
if (params.SpatialFlag)
    FreqStep = (FreqMax - FreqMin) / FreqPts;
else
    FreqStep = (FreqMax - FreqMin) / FreqPts / 2;
end

clear freqSweepStim
for t = 1:FreqPts
    freqSweepStim(t) = sin(2*pi*(FreqMin + t*FreqStep)*t*1e-4);
end

clear contrastSweepStim
for t = 1:ContrastPts
    contrastSweepStim(t) = 1 + sin(2*pi*t / ContrastPeriod) * (ContrastMin + (t / ContrastPts) * (ContrastMax-ContrastMin));
end
gapStim(1:InterPts) = 1;

clear IncDecStim;
IncDecStim(1:InterPts) = 1;
IncDecStim(InterPts+1:InterPts+StepPts) = 2;
IncDecStim(InterPts+StepPts+1:InterPts*2 + StepPts) = 1;
IncDecStim(2*InterPts+StepPts+1:2*InterPts+2*StepPts) = 0;
IncDecStim(2*InterPts+2*StepPts+1:3*InterPts+2*StepPts) = 1;

% analyze inc/dec ratio
if (params.SpatialFlag)
    incResp = MeanResponse(PrePts/2 + 1:PrePts + StepPts);
    decResp = MeanResponse(PrePts + InterPts + StepPts + 1 - PrePts/2:PrePts + InterPts + 2*StepPts);
else
    incResp = MeanResponse(3*InterPts/4 + 1:InterPts + StepPts);
    decResp = MeanResponse(2*InterPts + StepPts + 1 - InterPts/4:2*InterPts + 2*StepPts);
end
if (params.CellAttached)
    results.IncDecRatio = max(incResp) / min(decResp);
else
    results.IncDecRatio = -min(incResp) / max(decResp);
    results.MaxInc = -min(incResp);
    results.MaxDec = max(decResp);
    [maxVal, results.tPeakInc] = min(incResp);
    results.tPeakInc = (results.tPeakInc - InterPts/4) * 1e-4;
end

figure(3);
subplot(4, 1, 2);
tme = 1:length(incResp);
tme = (tme - InterPts/4) * 1e-4;
plot(tme, incResp,  'LineWidth', 2); hold on; plot(tme, decResp,  'LineWidth', 2);
ylabel('pA')
xlabel('sec')
axis tight;
if (params.SaveToIgor)
    CurFigPanel = gca;
    makeAxisStruct(CurFigPanel,'ChirpIncDec','test')
end

if (params.SpatialFlag)
    freqSweep = MeanResponse(PrePts + 2*InterPts + 2*StepPts + 1:PrePts + 2*InterPts + 2*StepPts + FreqPts);
    contrastSweep = MeanResponse(PrePts + 3*InterPts + 2*StepPts + FreqPts + 1:PrePts + 3*InterPts + 2*StepPts + FreqPts + ContrastPts);
else
    freqSweep = MeanResponse(3*InterPts + 2*StepPts + 1:3*InterPts + 2*StepPts + FreqPts);
    contrastSweep = MeanResponse(4*InterPts + 2*StepPts + FreqPts + 1:4*InterPts + 2*StepPts + FreqPts + ContrastPts);
end

figure(2); clf
subplot(2, 1, 1);
plot(freqSweep);
subplot(2, 1, 2);
plot(contrastSweep);
pause(1);

% analyze contrast sweep
clear Contrast ContrastResp
for cycle = 1:ContrastCycles
    Dat = contrastSweep((cycle-1)*ContrastPeriod+1:cycle*ContrastPeriod);
    ContrastResp(cycle) = -(min(Dat) - max(Dat));
    Contrast(cycle) = ContrastMin + (cycle-1)*(ContrastMax - ContrastMin) / ContrastCycles;
end

coef = [0.1 0.1];
fitcoef = nlinfitsome([false false], Contrast, ContrastResp ./ max(ContrastResp), @hill, coef);
fit = hill(fitcoef, Contrast);
figure(3); subplot(4, 1, 4);
plot(Contrast, ContrastResp ./ max(ContrastResp), 'o', Contrast, fit);
xlabel('contrast');
ylabel('norm response');
Indices = find(fit > 0.5);
if (Indices(1) == 1)
    Indices(1) = 2;
end
Indices2 = find(fit(Indices-1) < 0.5);
results.HalfMaxContrast = mean(Contrast(Indices(Indices2)));
if (params.SaveToIgor)
    CurFigPanel = gca;
    makeAxisStruct(CurFigPanel,'ChirpContrast','test')
end

if (params.CellAttached)
    results.IncDecRatioSine = max(contrastSweep) / min(contrastSweep);
else
    results.IncDecRatioSine = -min(contrastSweep) / max(contrastSweep);
end

% analyze frequency sweep
Indices = find(freqSweepStim > 0);
if (Indices(1) == 1);
    Indices(1) = 2;
end

Indices2 = find(freqSweepStim(Indices - 1) < 0);

ZeroCrossIndices = Indices(Indices2);

clear freq amp;
for indx = 1:length(ZeroCrossIndices)-1
    freq(indx) = FreqMin + ZeroCrossIndices(indx)*FreqStep;
    amp(indx) = max(freqSweep(ZeroCrossIndices(indx):ZeroCrossIndices(indx+1))) - min(freqSweep(ZeroCrossIndices(indx):ZeroCrossIndices(indx+1)));
end
coef = [1 -2];
fitcoef = nlinfitsome([false false], freq, amp ./ max(amp), @hill, coef);
fit = hill(fitcoef, freq);
figure(3); subplot(4, 1, 3);
plot(freq, amp/max(abs(amp)), 'o', freq, fit);
Indices = find(fit < 0.5);
Indices2 = find(fit(Indices-1) > 0.5);
if (length(Indices) >= length(Indices2))
    results.HalfAttenuatingFreq = mean(freq(Indices(Indices2)));
else
    results.HalfAttenuatingFreq = nan;        
end
xlabel('frequency (Hz)');
ylabel('norm response');
if (params.SaveToIgor)
    CurFigPanel = gca;
    makeAxisStruct(CurFigPanel,'ChirpFrequency','test')
end

if (params.SaveGraphs)
    fileName = strcat(strcat(titleString), '-');
    fileName = strcat(fileName, params.fileName);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
    saveas(gcf, strcat(fileName, '.pdf'));
end


