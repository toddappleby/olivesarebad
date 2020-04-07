%% initialize save structs
MNOFFParasol = struct();
MNONParasol = struct();
MNONSmooth = struct();
MNOFFSmooth = struct();

% MCS 

MCSOFFParasol = struct();
MCSONParasol = struct();
MCSONSmooth = struct();
MCSOFFSmooth = struct();
%% I'M DUMB
OFFParasolC = 0;
ONParasolC = 0;
ONSmoothC = 0;
OFFSmoothC = 0;

%% load FT first, give different name --- #1
cd('E:\Data Analysis_2020\2020_0406\')
expDate = dir;
load('20200406Bc1.mat')
frameTimings = epochs;
%% now load data -- #2

% load('20200114Bc2.mat')

uniqueProtocols = [];

for z = 1:size(epochs,2)
   list(z) = string(epochs(z).meta.displayName);
%    allEpochData(z,1:length(epochs(z).epoch)) = epochs(z).epoch;
end
%find 0s, make index, create new string
list
while ~isempty(list)
uniqueProtocols = [uniqueProtocols; list(1)];%OK<AGROW>
uniqueCheck = strcmp(uniqueProtocols(length(uniqueProtocols)),list);
newIndex = find(uniqueCheck==0);
list = list(newIndex);
% also store epochs this way
currentIndex = find(uniqueCheck==1);
% protocolData = allEpochData(currentIndex,:);
rawData{length(uniqueProtocols),1} = uniqueProtocols(length(uniqueProtocols));
% rawData{length(uniqueProtocols),2} = protocolData;
% allEpochData = allEpochData(newIndex,:);
end

uniqueProtocols = sort(uniqueProtocols) %#ok<NOPTS>

% SAVE PARAMS HERE:
cellType = 'OFF Parasol';
saveFlag = 0;
% cellNum = '1';
expNum= 1;
cellName = 'Bc1';
%% spike detection and save into protocol groups -- NOT IN USE

for s = 1:length(uniqueProtocols)
    currentProtocol = uniqueProtocols(s); 
    for b = 1:size(epochs,2)
        allEpochs(b,:) = epochs(b).epoch;
        
    end
end

%% frame timings and epoch organization -- #3
count = 0;
clear monitorStorage
clear epochStartTime
clear epochNumFrames
for f = 1:length(frameTimings)
 
     displayName = frameTimings(f).meta.displayName;
     
   
   if strcmp(displayName,'Motion And Noise')
       
    count = count+1;
    singleMonitorRun = frameTimings(f).epoch;
 monitorStorage(count,:) = singleMonitorRun;
 epochStartTime(count) = frameTimings(f).meta.epochTime;
 epochNumFrames(count) = frameTimings(f).meta.epochNum;
   end
end

[dates,DI] = sort(epochStartTime);  
monitorStorage = monitorStorage(DI,:);
epochNumFrames = epochNumFrames(DI);


%%  spike detection and epoch organization -- #4

 count = 0;
clear binnedStorage

clear epochStartTimeD

for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   
   if strcmp(displayName,'Motion And Noise')
        
        
count = count + 1;
        
sampleRate = 10000;
frameRate = 60;
binRate = 1000;
preTime = epochs(i).meta.preTime;
stimTime = epochs(i).meta.stimTime;
% stimTimes(count)=epochs(i).meta.stimTime;
stimOrig = stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
tailTime = epochs(i).meta.tailTime;
       
background = epochs(i).meta.backgroundClass;
      
prePts = preTime * 1e-3 * binRate;
stimPts = epochs(i).meta.stimTime * 1e-3 * binRate;
        
 singleEpoch = epochs(i).epoch;
%  binnedRow = binData(singleEpoch,binRate,sampleRate);
%   binnedRow = binData(singleEpoch(prePts+1:end),binRate,sampleRate);
% binnedStorage(count,:) = binnedRow(:)';
binnedStorage(count,:) = singleEpoch;
%  epochStorage(count,:) = singleEpoch;
 epochNumResponses(count) = epochs(i).meta.epochNum;
 
seedList(count) = epochs(i).meta.seed;

epochStartTimeD(count) = epochs(i).meta.epochTime;
% frameDwell(count) = epochs(i).meta.frameDwell;
noiseClass(count) = {epochs(i).meta.noiseClass};

        switch background
            case 'sequential'
                bgNum = 1;
            case 'random'
                bgNum = 2;
            case 'stationary'
                bgNum = 3;
        end
                bgClass(count,1) = bgNum;
     
       end

  end
     
spikeMatrix = zeros(size(binnedStorage,1),size(binnedStorage,2)/10);
spikeMatrixUnbinned = zeros(size(binnedStorage,1),size(binnedStorage,2));
psthMatrix = zeros(size(binnedStorage,1),size(binnedStorage,2)/10);

desiredSTD = 7; %arbitrary, works for most extracellular data

% binRate = 1000;
% response(o,:) = binSpikeCount(epoch(o,:)/sampleRate, binRate, sampleRate);
% response(o,:) = psth(response(o,:)*binRate,6+2/3,binRate,1);

%spike detection
for k = 1:size(binnedStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(binnedStorage(k,:), stimTime, ...
              desiredSTD);
          
          if isempty(spikes)
              disp('deleted epoch')
          else
               
        spikeMatrixUnbinned(k,:) = spikes;
%         psthMatrix2(k,:) = psth(spikeMatrix(k,:)*10000,6+2/3,10000,1);
        spikes = binSpikeCount(spikes/sampleRate, binRate, sampleRate);
        spikes = spikes';
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:)*binRate,6+2/3,binRate,1);
%         
%         spikeMatrix(k,:) = binSpikeCount(spikes/sampleRate, binRate, sampleRate);
%         psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
        
          end
end

%  analysisType = 'extracellular';
% for k = 1:size(binnedStorage,1)
%     
%    spikes = responseByType(epochStorage(k,:), analysisType, preTime, sampleRate);
%     spikeMatrix(k,:) = binSpikeCount(spikes/sampleRate, binRate, sampleRate);
%     psthMatrix(k,:) = psth(spikeMatrix(k,:)*binRate,6+2/3,binRate,1);
% 
% end

%NOTE: Below uses start time index to bring all data and metadata into
%register with frametime epochs

[datesD,DI2] = sort(epochStartTimeD);
spikeMatrix = spikeMatrix(DI2,:);
spikeMatrixUnbinned = spikeMatrixUnbinned(DI2,:);
epochNumResponses = epochNumResponses(DI2);
psthMatrix = psthMatrix(DI2,:);
 bgClass = bgClass(DI2);
 seedList = seedList(DI2);
 noiseClass = noiseClass(DI2);
 binaryIndex = find(contains(noiseClass,'binary'));
 gaussianIndex = find(contains(noiseClass,'gaussian'));
 bgClassOrig = bgClass;
 
%% reorganize so indices match between frametimes and responses (now done above..)

% [matches,frameI,stimI] = intersect(epochNumFrames,epochNumResponses);
% spikeMatrix = spikeMatrix(stimI,:);
% bgClass = bgClass(stimI);
% seedList = seedList(stimI);
% monitorStorage = monitorStorage(frameI,:);




%% Get stim frames -- #5 (standard, using Mike's frame time functions)
count = 0;
noiseFlag =1; %1 for gaussian
clear frameValues;
clear frameValuesAll;
clear response1;
%number of Frames, frame dwell, st dev, seed
% Note: I almost never change the noise contrast, but will need to deal
% with binary noise.  That can be done in the initial sort for the most
% part I think.
if noiseFlag == 1
noiseVars = struct();
noiseVars.type = 'gaussian';
noiseVars.contrast = 0.3333;

timings = [250,10000,250]; % AUTOMATE THIS LATER
frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage(gaussianIndex,:), 10000, binRate);
% frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage, 10000, 10000);
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,1,seedList(gaussianIndex));
else
    noiseVars = struct();
noiseVars.type = 'binary';
noiseVars.contrast = 0.3333;

timings = [250,10000,250]; % AUTOMATE THIS LATER
frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage(binaryIndex,:), 10000, binRate);
% frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage, 10000, 10000);
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,1,seedList(binaryIndex));
end
    
%% Max Turner frame timing stuff -- NOT IN USE
clear frameTimes
clear S
clear R
clear binnedStimulus
clear binnedResponse
stimFrames = round(frameRate * (stimOrig/1e3));
frameDwell = 1;
count = 0;
countSeq = 0;
linearFilterSeq=0;

% frameTimes = getFrameTiming(monitorStorage(t,:),2);
% testStim = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frameTimes,1,seedList);
frameTimesCell = cell(1,size(monitorStorage,1));

for t = 1:size(monitorStorage,1)
   
frameTimes = getFrameTiming(monitorStorage(t,:),2); %apparently 2 is for OLED
frameTimesCell{t} = frameTimes;
preFrames = frameRate*(preTime/1000);
firstStimFrameFlip = frameTimes(preFrames);
spikeSelect = spikeMatrixUnbinned(t,firstStimFrameFlip:end);

noiseStream = RandStream('mt19937ar', 'Seed', seedList(t));

spacerWidth = mean(diff(frameTimes));

    for tt = 1:floor(stimFrames/frameDwell)
        
        binnedStimulus(tt) = .3 * noiseStream.randn;
        binnedResponse(tt) = mean(spikeSelect(1,(round((tt-1)*spacerWidth + 1) : round(tt*spacerWidth))));
        
    end
S(t,:) = binnedStimulus;
R(t,:) = binnedResponse;
    
% mhtFrameValues(t,:) = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),binRate,mhtFrames,1,seedList(t));
end

testStim = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),10000,frameTimesCell,1,seedList);


SI=find(bgClass==1);
lfilterS = getLinearFilter(S(SI,:), R(SI,:), ...
'analysisType', 'revcorr', ...
'fourierCorrection', false, ...
'binRate', spacerWidth, ...
'filterTime', 0.5, ...
'frameRate', 60);

for h = 1:size(S,1)
    
    if bgClass(h) == 1
        countSeq = countSeq+1;
lfResp = spikeMatrixUnbinned(h,:);
lfStim = testStim(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterSeq = (linearFilterSeq*(countSeq-1) + lf)/countSeq;
sortNLSeq(countSeq) = h;
    end
end



% frameValuesAll = mhtFrameValues;
%% probably superfluous, but does take out 1 seeds
count = 0;
clear response1;
clear frameValues;
for g = 1:size(binnedStorage,1)
      if seedList(g) ~= 1
        
        count=count+1;
% frameValRow = getGaussianNoiseFrames((10000/1000) * 60,1,(1/3)*0.3,seedList(g));
% points per frame? 
% numP = round(binRate/60);
% frameValRow = ones(numP,1)*frameValRow(:)';
% frameValues(count,:) = frameValRow(:);

% response1(count,:) = spikeMatrixUnbinned(g,1:size(frameValuesAll,2));

% response1(count,:) = psthMatrix2(g,:);
frameValues(count,:) = frameValuesAll(g,:);
bgClass(count) = bgClass(g);

response1(count,:) = psthMatrix(g,:);
% frameValues(count,:) = testStim(g,:);
% bgClass(count) = bgClass(g);

      end
end
plotLngth = round(binRate*.5);
% response1(:,1:floor(binRate/2)) = 0;
% frameValues(:,1:floor(binRate/2)) = 0;
 bgClass(count+1:end)=[];
%% just use this for gaussian/binary - #6 (keeps 1 seeds for now)
plotLngth = round(binRate*.5);
clear response1
clear frameValues
bgClass = bgClassOrig;
if noiseFlag == 1
     response1 = psthMatrix(gaussianIndex,:);
 frameValues = frameValuesAll;
 bgClass = bgClass(gaussianIndex,:);
else
 response1 = psthMatrix(binaryIndex,:);
 frameValues = frameValuesAll;
 bgClass = bgClass(binaryIndex,:);
end
%% correlation stuff -- #7
clear seqParams
linearFilter = 0;
linearFilterSeq = 0;
linearFilterRandom = 0;
linearFilterStatic = 0;
linearFilterNew = 0;
countSeq =0;
countR = 0;
countStatic =0;
sortNLSeq=0;
sortNLRand=0;
sortNLStatic=0;
for h = 1:size(response1,1)
    
    if bgClass(h) == 1
        countSeq = countSeq+1;
lfResp = response1(h,:);
lfStim = frameValues(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterSeq = (linearFilterSeq*(countSeq-1) + lf)/countSeq;
sortNLSeq(countSeq) = h;


    elseif bgClass(h) == 2
        countR = countR+1;
lfResp = response1(h,:);
lfStim = frameValues(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterRandom = (linearFilterRandom*(countR-1) + lf)/countR;
sortNLRand(countR) = h;

    elseif bgClass(h) == 3
        countStatic = countStatic + 1;
lfResp = response1(h,:);
lfStim = frameValues(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterStatic = (linearFilterStatic*(countStatic-1) + lf)/countStatic;
sortNLStatic(countStatic) = h;

    end
        
end

%normalize to scale the filters 
linearFilterSeq = linearFilterSeq/norm(linearFilterSeq);
linearFilterRandom = linearFilterRandom/norm(linearFilterRandom);
linearFilterStatic = linearFilterStatic/norm(linearFilterStatic);

%fits for filters using Obsidian functions 
     seqParams = fitLinearFilterParams(linearFilterSeq,binRate);
     filterSeq = linearFilterFunction(seqParams,(1:plotLngth)/binRate);
     figure(7)
     plot((1:plotLngth)/binRate, filterSeq,'LineWidth',2,'Color','r')
     hold on
     
     randomParams = fitLinearFilterParams(linearFilterRandom,binRate);
     filterRand = linearFilterFunction(randomParams,(1:plotLngth)/binRate);
     plot((1:plotLngth)/binRate, filterRand,'LineWidth',2,'Color','b')
     
     staticParams = fitLinearFilterParams(linearFilterStatic,binRate);
     filterStatic = linearFilterFunction(staticParams,(1:plotLngth)/binRate);
     plot((1:plotLngth)/binRate, filterStatic,'LineWidth',2,'Color','k')
     set(gca,'xdir','reverse')
     legend('sequential','random','static','location','northwest')
     ylabel('weight')
     xlabel('time (s)')
     
%      figure(11)
%      plot((1:plotLngth)/binRate,linearFilterNew)
     

     
%% plotlf - raw filters #8

figure(6)
clf
plot((1:plotLngth)/binRate, linearFilterSeq(1:plotLngth),'Color','r')
hold on
line((1:plotLngth)/binRate, linearFilterRandom(1:plotLngth),'Color','b')
line((1:plotLngth)/binRate, linearFilterStatic(1:plotLngth),'Color','k')
set(gca,'xdir','reverse')
title('linear filters for each surround condition')
legend('Sequential','Random','Static','Location','northwest')
xlabel('time(s)')
ylabel('weight')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
export_fig 'noiseNL.pdf'

 if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
    OFFParasolC = OFFParasolC+1;
    
     MNOFFParasol.Exp(expNum).Cell(OFFParasolC).Name = cellName;
     MNOFFParasol.Exp(expNum).Cell(OFFParasolC).Date = epochs(1).meta.epochTime;
     MNOFFParasol.Exp(expNum).Cell(OFFParasolC).LF.X= ((1:plotLngth)/binRate);
     MNOFFParasol.Exp(expNum).Cell(OFFParasolC).LF.Seq = linearFilterSeq;
     MNOFFParasol.Exp(expNum).Cell(OFFParasolC).LF.Rand = linearFilterRandom;
     MNOFFParasol.Exp(expNum).Cell(OFFParasolC).LF.Static = linearFilterStatic;
 elseif saveFlag ==1 && strcmp(cellType,'ON Parasol')
    ONParasolC = ONParasolC+1;
    
    MNONParasol.Exp(expNum).Cell(ONParasolC).Name = cellName;
    MNONParasol.Exp(expNum).Cell(ONParasolC).Date = epochs(1).meta.epochTime;
    MNONParasol.Exp(expNum).Cell(ONParasolC).LF.X= (1:plotLngth)/binRate;
    MNONParasol.Exp(expNum).Cell(ONParasolC).LF.Seq = linearFilterSeq;
    MNONParasol.Exp(expNum).Cell(ONParasolC).LF.Rand = linearFilterRandom;
    MNONParasol.Exp(expNum).Cell(ONParasolC).LF.Static = linearFilterStatic;
 elseif saveFlag ==1 && strcmp(cellType,'ON Smooth')
    ONSmoothC = ONSmoothC+1;
     
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).Name = cellName;
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).Date = epochs(1).meta.epochTime;
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).LF.X= (1:plotLngth)/binRate;
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).LF.Seq = linearFilterSeq;
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).LF.Rand = linearFilterRandom;
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).LF.Static = linearFilterStatic;
 elseif saveFlag ==1 && strcmp(cellType,'OFF Smooth')
    OFFSmoothC = ONSmoothC+1;
     
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).Name = cellName;
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).Date = epochs(1).meta.epochTime;
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).LF.X= (1:plotLngth)/binRate;
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).LF.Seq = linearFilterSeq;
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).LF.Rand = linearFilterRandom;
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).LF.Static = linearFilterStatic;
 end
%% tests to see if Mike's filter function changes anything (doesn't seem to)
SI=find(bgClass==1);
lfilterS = getLinearFilter(frameValues(SI,:), response1(SI,:), ...
'analysisType', 'revcorr', ...
'fourierCorrection', false, ...
'binRate', binRate, ...
'filterTime', 0.5, ...
'frameRate', 60);

RI = find(bgClass==2);

lfilterR = getLinearFilter(frameValues(RI,:), response1(RI,:), ...
'analysisType', 'revcorr', ...
'fourierCorrection', false, ...
'binRate', binRate, ...
'filterTime', 0.5, ...
'frameRate', 60);

StI = find(bgClass==3);

lfilterSt = getLinearFilter(frameValues(StI,:), response1(StI,:), ...
'analysisType', 'revcorr', ...
'fourierCorrection', false, ...
'binRate', binRate, ...
'filterTime', 0.5, ...
'frameRate', 60);

figure(12)
plot(lfilterS(1:5000),'Color','r')
hold on
plot(lfilterR(1:5000),'Color','b')
plot(lfilterSt(1:5000),'Color','k')
set(gca,'xdir','reverse')
%% start NL 
stimExtent = prePts + (1:stimPts);
xaxis=0;
clear respStore;
clear pred2Store;

nonlinearityBins = 100;

%Convolve stimulus with filter to get generator signal.
seqInd = find(bgClass==1); seqFV = frameValues(sortNLSeq,:);
seqR = response1(sortNLSeq,:);

 for p = 1:size(seqFV,1)

yaxis=0;
xaxis = 0;
pred=0;
resp=0;


pred = ifft(fft(seqFV(p,:)) .* fft(linearFilterSeq(:)'));
resp = binData(seqR(p,:), 60, binRate);
yaxis = [yaxis, resp(:)'];
yaxisStore(p,:) = yaxis;
pred2Store(p,:)=pred;
pred = binData(pred,60,binRate);
pred=pred(:)';
predStore(p,:)=pred;
respStore(p,:)=[seqR(p,:) zeros(1,100)];

xaxis = [xaxis, pred(1 : length(resp))];
[a, b] = sort(xaxis(:));
 xSort = a;
 ySort = yaxis(b);
   % Bin the data.
 valsPerBin = floor(length(xSort) / nonlinearityBins);
 xBin = mean(reshape(xSort(1 : nonlinearityBins*valsPerBin),valsPerBin,nonlinearityBins));
 yBin = mean(reshape(ySort(1 : nonlinearityBins*valsPerBin),valsPerBin,nonlinearityBins));
 
 xBinHolder(p,:) = xBin;
 yBinHolder(p,:) = yBin;
 
 end

 huhx = mean(xBinHolder);
 huhy = mean(yBinHolder);
 figure(8);
 clf
 
 plot(huhx*10,huhy*10000,'Color','r')
% axis([-.5 1 0 100])
 hold on
%  predStore(:,stimExtent(501:end))
% [xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(1:9770)),respStore(:,stimExtent(1:9770)),nonlinearityBins);
[xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(:)),respStore(:,stimExtent(:)),nonlinearityBins);
 nlParams = fitNonlinearityParams(xfBin, yfBin);
 nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),1000);
 figure(9)
 plot(nlX,outputNonlinearity(nlParams,nlX),'Color','r','LineWidth',2)
 hold on
 
 if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
    MNOFFParasol.Exp(expNum).Cell(OFFParasolC).NL.X = nlX;
    MNOFFParasol.Exp(expNum).Cell(OFFParasolC).NL.Seq = outputNonlinearity(nlParams,nlX);
 elseif saveFlag ==1 && strcmp(cellType,'ON Parasol')
    MNONParasol.Exp(expNum).Cell(ONParasolC).NL.X = nlX;
    MNONParasol.Exp(expNum).Cell(ONParasolC).NL.Seq = outputNonlinearity(nlParams,nlX);
 elseif saveFlag ==1 && strcmp(cellType,'ON Smooth')
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).NL.X = nlX;
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).NL.Seq = outputNonlinearity(nlParams,nlX);
 elseif saveFlag ==1 && strcmp(cellType,'OFF Smooth')
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).NL.X = nlX;
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).NL.Seq = outputNonlinearity(nlParams,nlX);
 end
 
%% random NL 

xaxis=0;

nonlinearityBins = 100;

%Convolve stimulus with filter to get generator signal.
randInd = find(bgClass==2); randFV = frameValues(sortNLRand,:);
randR = response1(sortNLRand,:);

 for p = 1:size(randFV,1)

yaxis=0;
xaxis = 0;
pred=0;
resp=0;
pred = ifft(fft(randFV(p,:)) .* fft(linearFilterRandom(:)'));
resp = binData(randR(p,:), 60, binRate);
yaxis = [yaxis, resp(:)'];
yaxisStore(p,:) = yaxis;
pred2Store(p,:)=pred;
pred = binData(pred,60,binRate);
pred=pred(:)';
predStore(p,:)=pred;
respStore(p,:)=[randR(p,:) zeros(1,100)];

xaxis = [xaxis, pred(1 : length(resp))];
[a, b] = sort(xaxis(:));
 xSort = a;
 ySort = yaxis(b);
   % Bin the data.
 valsPerBin = floor(length(xSort) / nonlinearityBins);
 xBin = mean(reshape(xSort(1 : nonlinearityBins*valsPerBin),valsPerBin,nonlinearityBins));
 yBin = mean(reshape(ySort(1 : nonlinearityBins*valsPerBin),valsPerBin,nonlinearityBins));
 
 xBinHolder(p,:) = xBin;
 yBinHolder(p,:) = yBin;
 
 end

 huhx = mean(xBinHolder);
 huhy = mean(yBinHolder);
 %axis([0 600 0 10*10^-3])
 figure(8)
 plot(huhx*10,huhy*10000,'Color','b')
 [xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(:)),respStore(:,stimExtent(:)),nonlinearityBins);
%  [xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(1:9770)),respStore(:,stimExtent(1:9770)),nonlinearityBins);
 nlParams = fitNonlinearityParams(xfBin, yfBin);
 nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),1000);
 figure(9)
 plot(nlX,outputNonlinearity(nlParams,nlX),'Color','b','LineWidth',2)
 hold on
 
 if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
   MNOFFParasol.Exp(expNum).Cell(OFFParasolC).NL.Rand = outputNonlinearity(nlParams,nlX);
   elseif saveFlag ==1 && strcmp(cellType,'ON Parasol')
   MNONParasol.Exp(expNum).Cell(ONParasolC).NL.Rand = outputNonlinearity(nlParams,nlX);
   elseif saveFlag ==1 && strcmp(cellType,'ON Smooth')
   
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).NL.Rand = outputNonlinearity(nlParams,nlX);
   elseif saveFlag ==1 && strcmp(cellType,'OFF Smooth')
    
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).NL.Rand = outputNonlinearity(nlParams,nlX);
 end
%% stationary NL 

xaxis=0;

nonlinearityBins = 100;

%Convolve stimulus with filter to get generator signal.
staticInd = find(bgClass==3); staticFV = frameValues(sortNLStatic,:);
staticR = response1(sortNLStatic,:);

 for p = 1:size(staticFV,1)

yaxis=0;
xaxis = 0;
pred=0;
resp=0;
pred = ifft(fft(staticFV(p,:)) .* fft(linearFilterStatic(:)'));
resp = binData(staticR(p,:), 60, binRate);
yaxis = [yaxis, resp(:)'];
yaxisStore(p,:) = yaxis;
pred2Store(p,:)=pred;
pred = binData(pred,60,binRate);
pred=pred(:)';
predStore(p,:)=pred;
respStore(p,:)=[staticR(p,:) zeros(1,100)];

xaxis = [xaxis, pred(1 : length(resp))];
[a, b] = sort(xaxis(:));
 xSort = a;
 ySort = yaxis(b);
   % Bin the data.
 valsPerBin = floor(length(xSort) / nonlinearityBins);
 xBin = mean(reshape(xSort(1 : nonlinearityBins*valsPerBin),valsPerBin,nonlinearityBins));
 yBin = mean(reshape(ySort(1 : nonlinearityBins*valsPerBin),valsPerBin,nonlinearityBins));
 
 xBinHolder(p,:) = xBin;
 yBinHolder(p,:) = yBin;
 
 end

 huhx = mean(xBinHolder);
 huhy = mean(yBinHolder);
 figure(8)
 plot(huhx*10,huhy*10000,'Color','k')
%   axis([0 400 0 4*10^-3])
 hold on
 legend('sequential','random','static','location','northwest')
 title('NL for each surround condition')
 xlabel('input','FontSize',18)
 ylabel('spikes/s','FontSize',18)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
export_fig 'noiseNL.pdf' -append

[xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(:)),respStore(:,stimExtent(:)),nonlinearityBins);
%  [xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(1:9770)),respStore(:,stimExtent(1:9770)),nonlinearityBins);
 nlParams = fitNonlinearityParams(xfBin, yfBin);
 nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),1000);
 figure(9)
 plot(nlX,outputNonlinearity(nlParams,nlX),'Color','k','LineWidth',2)
  legend('sequential','random','static','location','northwest')
  xlabel('input')
  ylabel('spikes/s')

if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
    MNOFFParasol.Exp(expNum).Cell(OFFParasolC).NL.Static = outputNonlinearity(nlParams,nlX);
elseif saveFlag == 1 && strcmp(cellType,'ON Parasol')
    MNONParasol.Exp(expNum).Cell(ONParasolC).NL.Static = outputNonlinearity(nlParams,nlX);
elseif saveFlag ==1 && strcmp(cellType,'ON Smooth')
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).NL.Static = outputNonlinearity(nlParams,nlX);
elseif saveFlag ==1 && strcmp(cellType,'OFF Smooth')
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).NL.Static = outputNonlinearity(nlParams,nlX);
end

%  x=-300:400;
% nlfun = @(p,x)(p(1)*normcdf(p(2)*x+p(3),0,1));
% p = nlinfit(huhx,huhy,nlfun,[5 50 0]);
% 
% 
% 
% plot(xBin,yBin,'.');
% plot(x,nlfun(p,x));
% hold off;

%% final save
cd('/Users/toddappleby/Documents/Data/Clarinet Exports/SavedData')
save('MotionandNoise.mat','MNOFFParasol')
save('MotionandNoise.mat','-append','MNONParasol')
save('MotionandNoise.mat',',-append','MNONSmooth')
save('MotionandNoise.mat','-append','MNOFFSmooth')
%% Contrast Response Function
count = 0;
clear epochStorage
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   egLabel = epochs(i).meta.epochGroupLabel;
%    analysisType = epochs(i).meta.onlineAnalysis;
   if strcmp(displayName,'Contrast Response Spot') && ~strcmp(recordingTechnique,'whole-cell') && strcmp(egLabel,'Control')
        
        
count = count + 1;
epochStorage(count,:) = epochs(i).epoch;
contrast(count) = epochs(i).meta.contrast;
preTime = epochs(i).meta.preTime;
stimTime = epochs(i).meta.stimTime;
tailTime = epochs(i).meta.tailTime;
   end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 7;
spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end

% stimTime = epochs(1).meta.stimTime;
% tailTime = epochs(1).meta.tailTime;
timings = [preTime stimOrig tailTime];
saveStuff = struct();
saveStuff.flag = saveFlag;
saveStuff.cellType = cellType;
saveStuff.expNum = expNum;
saveStuff.OFFParasolC = OFFParasolC;
saveStuff.ONParasolC = ONParasolC;
if saveFlag ==1
[MNOFFParasol,MNONParasol] = CRF(spikeMatrix,psthMatrix,contrast,timings,saveStuff,MNOFFParasol,MNONParasol);
else
    CRF(spikeMatrix,psthMatrix,contrast,timings,saveStuff,[],[])
end


%% MTF spots and annuli 

splitFactors = ["stimulusClass","temporalFrequency","temporalClass","radius"];
runSplitter = ["radius"]; 

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage 

stringComparer = [];
stringComparer = string(stringComparer);
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;

   if strcmp(displayName,'S MT Fspot') && strcmp(egLabel,'Control')
        
        count = count + 1;
        for s = 1:length(splitFactors)
            if ischar(getfield(epochs(i).meta,splitFactors(s)))
                splitCell{2,s} = [splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];    
            else
                splitCell{2,s}= [splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
            end
        end
        
        count
        
%         width(count) = epochs(i).meta.barWidth;
tfreq(count) = epochs(i).meta.temporalFrequency;
        epochStorage(count,:) = epochs(i).epoch;
%         apertureRadius(count) = epochs(i).meta.apertureRadius;
%         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
%         angleO(count) = epochs(i).meta.orientation;
   end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 5;
spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));


% replacer = zeros(size(splitCell{2,1));

%DON'T GIVE A FUCK ABOUT A STRING!!!!

for m = 1:size(splitCell,2)-1
    if strcmp(class(splitCell{2,m}),'string')
        uniqueStrings = unique(splitCell{2,m});
        replacer = zeros(size(splitCell{2,m}));
        
        for n = 1:length(uniqueStrings)
           strInd = find(splitCell{2,m}==uniqueStrings(n));
           replacer(strInd) = n;
        end
        holdList(:,m) = splitCell{2,m};
        splitCell{2,m}=replacer;
    end
end
    

%SORT THIS SHIT BABY!!!!!!
allSets=zeros(size(splitCell{2,1},2),size(splitCell,2)-1);

for o = 1:size(splitCell,2)-1
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end




combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;

% while size(comboOut,1) ~= size(combos,1)
while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(istherenobetterway,2));
    combos(1,:) = [];
end


for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end
timings = [preTime stimOrig tailTime];
params = struct();
params.saveGraph = 0;
params.stimName = 'mTF';
params.tfreq = unique(tfreq);
% saveGraph = 0;
% stimName = 'Grating';
mTF(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);

%% Bar centering 

countx = 0;
county = 0;
clear epochStorage;
epochStorageX = [];
epochStorageY = [];
YspikeMatrix = [];
XspikeMatrix = [];
positionX = [];
positionY = [];

for i = 1:length(epochs)
   displayName = epochs(i).meta.displayName;
if strcmp(displayName,'Bar Centering')
   
       searchAxis = epochs(i).meta.searchAxis;
   if strcmp(searchAxis,'xaxis')
       
        countx = countx + 1;
        epochStorageX(countx,:) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
        positionX(countx) = epochs(i).meta.position;
        tFrequency = epochs(i).meta.temporalFrequency;
        
   elseif strcmp(searchAxis,'yaxis')
       
        county = county + 1;
        epochStorageY(county,:) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
        positionY(county) = epochs(i).meta.position;
        tFrequency = epochs(i).meta.temporalFrequency;
        
   end
end
   
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 5;
if ~isempty(epochStorageX)
    spikeMatrixX = zeros(size(epochStorageX,1),size(epochStorageX,2));
    psthMatrixX = zeros(size(epochStorageX,1),size(epochStorageX,2));

     disp('x') 
    for k = 1:size(epochStorageX,1)

        positionX(k)
    [spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorageX(k,:), stimTime, ...
                  desiredSTD);
              if isempty(spikes)
                  disp('deleted epoch')
              else
            XspikeMatrix(k,:) = spikes;
            XpsthMatrix(k,:) = psth(XspikeMatrix(k,:),6+2/3,sampleRate,1);
              end
    end
end

if ~isempty(epochStorageY)
    spikeMatrixY = zeros(size(epochStorageY,1),size(epochStorageY,2));
    psthMatrixY = zeros(size(epochStorageY,1),size(epochStorageY,2));

    disp('y')
    for k = 1:size(epochStorageY,1)
        positionY(k)

    [spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorageY(k,:), stimTime, ...
                  desiredSTD);
              if isempty(spikes)
                  disp('deleted epoch')
              else
            YspikeMatrix(k,:) = spikes;
            YpsthMatrix(k,:) = psth(YspikeMatrix(k,:),6+2/3,sampleRate,1);
              end
    end
end

% stimTime = epochs(1).meta.stimTime;
% tailTime = epochs(1).meta.tailTime;
timings = [preTime stimOrig tailTime];

centeringBars(XspikeMatrix,YspikeMatrix,positionX,positionY,timings,tFrequency,sampleRate)
%% Motion center surround
 

count = 0;
clear epochStorage 
clear bgClass
clear centerClass
clear contrastState

for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   egLabel = epochs(i).meta.epochGroupLabel;
   if strcmp(displayName,'Motion Center Surround') && ~strcmp(recordingTechnique,'whole-cell') && strcmp(egLabel,'Control')
        width = epochs(i).meta.surroundBarWidth;
        if width <= 100
        protocolExample = i;
        count = count + 1;
        epochStorage(count,:) = epochs(i).epoch;
        background = epochs(i).meta.backgroundClass;
        contrast = epochs(i).meta.contrast;
        center = epochs(i).meta.centerClass;
        switch background
            case 'sequential'
                bgNum = 1;
            case 'random'
                bgNum = 2;
            case 'stationary'
                bgNum = 3;
        end
        switch center
            case 'sequential'
                cNum = 1;
            case 'random'
                cNum = 2;
            case 'sequential180'
                cNum = 3;
        end
        bgClass(count,1) = bgNum;
        centerClass(count,1) = cNum; 
        contrastState(count,1) = contrast;
        end
   end
end

saveGraph = 0;
[meanCollection, meanKey] = motionCS(epochs,bgClass,centerClass,contrastState,epochStorage,protocolExample,saveGraph);



if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
    OFFParasolC = OFFParasolC + 1;
CMSOFFParasol(OFFParasolC,:) = meanCollection;
    elseif saveFlag == 1 && strcmp(cellType,'ON Parasol')
            ONParasolC = ONParasolC + 1;
            CMSONParasol(ONParasolC,:) = meanCollection;
        elseif saveFlag == 1 && strcmp(cellType,'ON Smooth')
                ONSmoothC = ONSmoothC + 1;
                CMSONSmooth(ONSmoothC,:) = meanCollection;
            elseif saveFlag == 1 && strcmp(cellType,'OFF Smooth')
                    OFFSmoothC = OFFSmoothC + 1;
                    CMSOFFSmooth(OFFSmoothC,:) = meanCollection;
end


if saveFlag == 1
testSave = strcat(expDate(3).name(1:8),'_CMS');
    cd('/Users/toddappleby/Documents/Data/Clarinet Exports/SavedData')
    save(testSave,'CMSOFFParasol')
    save(testSave,'-append','CMSONParasol')
    save(testSave,'-append','CMSONSmooth')
    save(testSave,'-append','CMSOFFSmooth')
end


%% LED Pulse

%this isn't as flexible because character vectors are hellworld

splitFactors = ["led"];
%NOTE: ISOMERIZATIONS PLEASE 

splitCell = cell(2,length(splitFactors));

for g = 1:size(splitCell,2)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage;

for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   if strcmp(displayName,'Led Pulse') && ~strcmp(recordingTechnique,'whole-cell') && ~strcmp(recordingTechnique,'EXCITATION') && ~strcmp(recordingTechnique,'INHIBITION')
        
        for s = 1:length(splitFactors)
            
        splitCell{2,s}=[splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];
        end
       
        count = count + 1;
%         intensity(count) = epochs(i).meta.intensity;

        epochStorage(count,:) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
      
   end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 4.5;
spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

allSets=zeros(size(splitCell{2,1},2),size(splitCell,2)-1);
clear allSets
for o = 1:size(splitCell,2)
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end

combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;

% while size(comboOut,1) ~= size(combos,1)
while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(comboOut,2));
    combos(1,:) = [];
end


for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch') 
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end

params = struct();
params.saveGraph =0;
params.stimName = 'LED';
timings = [preTime stimOrig tailTime];
LED(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% Expanding Spots

splitFactors = ["spotIntensity","backgroundIntensity","currentSpotSize"];
runSplitter = ["currentSpotSize"]; %might just set this manually because not gonna change? number could change if the protocol run was cut short
%this is why Greg did this the way he did .... maybe.  gonna have to figure
%this out by TIME??? which seems rough

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage 


for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
   if strcmp(displayName,'Expanding Spots') && strcmp(egLabel,'Control')
        
        count = count + 1;
        for s = 1:length(splitFactors)
        splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
        end
        
%         width(count) = epochs(i).meta.barWidth;
        epochStorage(count,:) = epochs(i).epoch;
%         apertureRadius(count) = epochs(i).meta.apertureRadius;
%         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
%         angleO(count) = epochs(i).meta.orientation;
   end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 4;
spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

%SORT THIS SHIT BABY!!!!!!
allSets=zeros(size(splitCell{2,1},2),size(splitCell,2)-1);

for o = 1:size(splitCell,2)-1
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end


combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;

% while size(comboOut,1) ~= size(combos,1)
while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(istherenobetterway,2));
    combos(1,:) = [];
end


for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end
timings = [preTime stimOrig tailTime];
params = struct();
params.saveGraph = 0;
params.stimName = 'expanding';
% saveGraph = 0;
% stimName = 'Grating';
simpleSpots(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);

%% Drifting Grating

splitFactors = ["apertureRadius","barWidth","temporalFrequency","orientation"];
runSplitter = ["orientations"]; %might just set this manually because not gonna change? number could change if the protocol run was cut short
%this is why Greg did this the way he did .... maybe.  gonna have to figure
%this out by TIME??? which seems rough

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage 
clear width
clear apertureRadius
clear temporalFrequency
clear angle0

for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
   if strcmp(displayName,'Grating DSOS') && strcmp(egLabel,'Control')
        
        count = count + 1;
        for s = 1:length(splitFactors)
        splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
        end
        
%         width(count) = epochs(i).meta.barWidth;
        epochStorage(count,:) = epochs(i).epoch;
%         apertureRadius(count) = epochs(i).meta.apertureRadius;
%         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
%         angleO(count) = epochs(i).meta.orientation;
   end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 5;
spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

%SORT THIS SHIT BABY!!!!!!
allSets=zeros(size(splitCell{2,1},2),size(splitCell,2)-1);

for o = 1:size(splitCell,2)-1
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end


combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;

% while size(comboOut,1) ~= size(combos,1)
while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==3);
    combos(1,:) = [];
end


for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end
timings = [preTime stimOrig tailTime];
params = struct();
params.saveGraph = 0;
params.stimName = 'Grating';
% saveGraph = 0;
% stimName = 'Grating';
orientedStim(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% Oriented Bars

splitFactors = ["intensity","backgroundIntensity","barSize","orientation"];
%NOTE: last split is always X axis.  I think this is helpful because can be
%specified by length function and don't need to cary another thing to
%processing function

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end


count = 0;
clear epochStorage;

for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   if strcmp(displayName,'Oriented Bars') && ~strcmp(recordingTechnique,'whole-cell') && ~strcmp(recordingTechnique,'EXCITATION') && ~strcmp(recordingTechnique,'INHIBITION')
        
        for s = 1:length(splitFactors)
        splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
        end
       
        count = count + 1;
        intensity(count) = epochs(i).meta.intensity;
        epochStorage(count,:) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
      
   end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 4.5;
spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

allSets=zeros(size(splitCell{2,1},2),size(splitCell,2)-1);

for o = 1:size(splitCell,2)-1
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end

combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;

% while size(comboOut,1) ~= size(combos,1)
while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(comboOut,2));
    combos(1,:) = [];
end


for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch') 
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end

params = struct();
params.saveGraph =0;
params.stimName = 'bars';
timings = [preTime stimOrig tailTime];
orientedStim(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% MOVING BAR


splitFactors = ["intensity","backgroundIntensity","barSize","orientation"];
%NOTE: last split is always X axis.  I think this is helpful because can be
%specified by length function and don't need to cary another thing to
%processing function

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage;

for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   
   if strcmp(displayName,'Moving Bar') && ~strcmp(recordingTechnique,'whole-cell') && ~strcmp(recordingTechnique,'EXCITATION') && ~strcmp(recordingTechnique,'INHIBITION') 
        
        for s = 1:length(splitFactors)
        splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
        end
       
        count = count + 1;
        intensity(count) = epochs(i).meta.intensity;
        epochStorage(count,:) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
      
   end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 6;
spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

allSets=zeros(size(splitCell{2,1},2),size(splitCell,2)-1);

for o = 1:size(splitCell,2)-1
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end

combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;

% while size(comboOut,1) ~= size(combos,1)
while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(comboOut,2));
    combos(1,:) = [];
end


for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch') 
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end

params = struct();
params.saveGraph =0;
%grating,bars
params.stimName = 'moving bar';
timings = [preTime stimOrig tailTime];
orientedStim(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);

%% processing functions

function LED(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

sampleRate = 10000;
stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;

figure(11)
subplot(2,1,1)

xvals = linspace(0,stimOff/10,7500);

plot(xvals,mean(psthMatrix(indexHolder{2,1},:)))
xline(stimStart/10,'LineStyle','--','Color','r','LineWidth',2)
xline(stimEnd/10,'LineStyle','--','Color','r','LineWidth',2)
ylabel(indexHolder{1,1})
xlabel('time (ms)')
subplot(2,1,2)
plot(xvals,mean(psthMatrix(indexHolder{2,2},:)))
xline(stimStart/10,'LineStyle','--','Color','m','LineWidth',2)
xline(stimEnd/10,'LineStyle','--','Color','m','LineWidth',2)
xlabel('time (ms)')
ylabel(indexHolder{1,2})


end
    
function [meanCollection, meanKey] = motionCS(epochs,bgClass,centerClass,contrastState,epochStorage,protocolExample,saveGraph)

sampleRate = 10000;
binRate = 1000;
preTime = epochs(protocolExample).meta.preTime;
stimTime = epochs(protocolExample).meta.stimTime;
stimOrig = stimTime;
motionTime = epochs(protocolExample).meta.motionTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 7;
% spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
% psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

for k = 1:size(epochStorage,1)
    bgClass(k,1);
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else
              
        spikes = binSpikeCount(spikes, binRate, sampleRate);
        spikes = spikes';
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:)*binRate,6+2/3,binRate,1);
%         spikeMatrix(k,:) = spikes;
%         psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end

size(spikeMatrix)
size(psthMatrix)


% sequenceOut = mean(psthMatrix(seq,:))';
% randomOut = mean(psthMatrix(rand,:))';
% staticOut = mean(psthMatrix(static,:))';
% xAxis = linspace(0,2000,20000);
% xAxis = xAxis';

% outTable = table(sequenceOut, randomOut, staticOut,xAxis);
% writetable(outTable,'outTable.xlsx','sheet',1);

% figure(26)
% clf
% plot(mean(psthMatrix(seq,1.16*10^4:1.28*10^4)),'Color','r')
% hold on
% plot(mean(psthMatrix(rand,1.16*10^4:1.28*10^4)),'Color','b')
% plot(mean(psthMatrix(static,1.16*10^4:1.28*10^4)),'Color','k')

% 
% spikeSeq = sum(spikeMatrix(seq,1.15*10^4:1.3*10^4),2);
% seqAvg = mean(spikeSeq)/.15;
% spikeRand = sum(spikeMatrix(rand,1.15*10^4:1.3*10^4),2);
% randAvg = mean(spikeRand)/.15;
% spikeStatic = sum(spikeMatrix(static,1.15*10^3:1.3*10^4),2);
% staticAvg = mean(spikeStatic)/.15;
% 
% maxSeq = max(mean(psthMatrix(seq,1.16*10^4:1.28*10^4)));
% maxRand = max(mean(psthMatrix(rand,1.16*10^4:1.28*10^4)));
% maxStatic = max(mean(psthMatrix(static,1.16*10^4:1.28*10^4)));


% spikeSeq = sum(psthMatrix(seq,.25*10^4:1*10^4),2);
% seqAvg = mean(spikeSeq)/.75;
% spikeRand = sum(spikeMatrix(rand,.25*10^4:1*10^4),2);
% randAvg = mean(spikeRand)/.75;
% spikeStatic = sum(spikeMatrix(static,.25*10^3:1*10^4),2);
% staticAvg = mean(spikeStatic)/.75;

% figure
% clf
%  parasolY=[seqAvg randAvg staticAvg];
%  smoothY = [seqAvg randAvg staticAvg];
%    dsY = [smoothY;0 0 0];
%   bar(dsY) 
% maxS = [maxSeq maxRand maxStatic];
% bar(dsY)

% 
seq = find(bgClass==1);
rand = find(bgClass==2);
static = find(bgClass==3);
 
posContrast = find(contrastState==1);
 negContrast = find(contrastState==-1);

posSeq = find(bgClass(posContrast) == 1);
posRand = find(bgClass(posContrast) == 2);
posStatic = find(bgClass(posContrast) == 3);

negSeq = find(bgClass(negContrast) == 1);
negRand = find(bgClass(negContrast) == 2);
negStatic = find(bgClass(negContrast) == 3);

centerSeq = find(centerClass((posContrast(posSeq)))==1);
centerRand = find(centerClass((posContrast(posSeq)))==2);
center180 = find(centerClass((posContrast(posSeq)))==3);

centerRSeq = find(centerClass((posContrast(posRand)))==1);
centerRRand = find(centerClass((posContrast(posRand)))==2);
centerR180 = find(centerClass((posContrast(posRand)))==3);

centerSeqNeg = find(centerClass((negContrast(negSeq)))==1);
centerRandNeg = find(centerClass((negContrast(negSeq)))==2);
center180Neg = find(centerClass((negContrast(negSeq)))==3);

centerRSeqNeg = find(centerClass((negContrast(negRand)))==1);
centerRRandNeg = find(centerClass((negContrast(negRand)))==2);
centerR180Neg = find(centerClass((negContrast(negRand)))==3);


returnStruct = struct();
returnStruct.surround.PSeq = mean(psthMatrix(posContrast(posSeq),:));
returnStruct.surround.PRand = mean(psthMatrix(posContrast(posRand),:));
returnStruct.surround.NSeq = mean(psthMatrix(negContrast(negSeq),:));
returnStruct.surround.NRand = mean(psthMatrix(negContrast(negRand),:));
returnStruct.center.PCS = mean(psthMatrix(posContrast(centerSeq),:));
returnStruct.center.PCR = mean(psthMatrix(posContrast(centerRand),:));
returnStruct.center.PCSt = mean(psthMatrix(posContrast(center180),:));
returnStruct.center.rPCS = mean(psthMatrix(posContrast(centerRSeq),:));
returnStruct.center.rPCR = mean(psthMatrix(posContrast(centerRRand),:));
returnStruct.center.rPCSt = mean(psthMatrix(posContrast(centerR180),:));
returnStruct.center.NCS = mean(psthMatrix(negContrast(centerSeqNeg),:));
returnStruct.center.NCR = mean(psthMatrix(negContrast(centerRandNeg),:));
returnStruct.center.NCSt = mean(psthMatrix(negContrast(center180Neg),:));
returnStruct.center.rNCS = mean(psthMatrix(negContrast(centerRSeqNeg),:));
returnStruct.center.rNCR = mean(psthMatrix(negContrast(centerRRandNeg),:));
returnStruct.center.rNCSt = mean(psthMatrix(negContrast(centerR180Neg),:));

motionStart = 100; % use these to extend/shorten spike count window
motionEnd = 150;

preTime = preTime/10;
motionTime=motionTime/10;

% figure(20)
% plot(mean(psthMatrix(negContrast(negSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd)))),'Color','r')
% hold on
% plot(mean(psthMatrix(negContrast(negRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd)))),'Color','b')

figure(10)
title('means')



% 
% a=sum(psthMatrix(posContrast(posSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% b=sum(psthMatrix(posContrast(posRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% 
% PCS = mean(a)/.25;
% PCSpeak = max(a);
% PCSerror = sem(a);
% 
% PCR = mean(b)/.25;
% PCRpeak = max(b);
% PCRerror = sem(b);
% 
% c=sum(psthMatrix(negContrast(negSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% d=sum(psthMatrix(negContrast(negRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% 
% NCS = mean(c)/.25;
% NCSpeak = max(c);
% NCSerror=sem(c);
% 
% NCR = mean(d)/.25;
% NCRpeak = max(d);
% NCRerror = sem(d);
% 
% e=sum(psthMatrix(posContrast(centerSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% f=sum(psthMatrix(posContrast(centerRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% g=sum(psthMatrix(posContrast(center180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% 
% CPCS= mean(e)/.25;
% CPCSpeak = max(e);
% CPCSerror = sem(e);
% CPCR = mean(f)/.25;
% CPCRpeak = max(f);
% CPCRerror = sem(f);
% CPC180 = mean(g)/.25;
% CPC180peak = max(g);
% CPC180error = sem(g);
% 
% h=sum(psthMatrix(posContrast(centerRSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% i=sum(psthMatrix(posContrast(centerRRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% j=sum(psthMatrix(posContrast(centerR180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% 
% rCPCS= mean(h)/.25;
% rCPCSpeak = max(h);
% rCPCSerror = sem(h);
% rCPCR = mean(i)/.25;
% rCPCRpeak = max(i);
% rCPCRerror = sem(i);
% rCPC180 = mean(j)/.25;
% rCPC180peak = max(j);
% rCPC180error = sem(j);
% 
% k=sum(psthMatrix(negContrast(centerSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% l=sum(psthMatrix(negContrast(centerRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% m=sum(psthMatrix(negContrast(center180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% 
% CNCS= mean(k)/.25;
% CNCSpeak = max(k);
% CNCSerror = sem(k);
% CNCR = mean(l)/.25;
% CNCRpeak = max(l);
% CNCRerror = sem(l);
% CNC180 = mean(m)/.25;
% CNC180peak = max(m);
% CNC180error = sem(m);
% 
% n=sum(psthMatrix(negContrast(centerRSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% o=sum(psthMatrix(negContrast(centerRRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% p=sum(psthMatrix(negContrast(centerR180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
% 
% rCNCS= mean(n)/.25;
% rCNCSpeak = max(n);
% rCNCSerror = sem(n);
% rCNCR = mean(o)/.25;
% rCNCRpeak = max(o);
% rCNCRerror = sem(o);
% rCNC180 = mean(p)/.25;
% rCNC180peak=max(p);
% rCNC180error = sem(p);


preTime = preTime/10;
motionTime=motionTime/10;

% motionStart = 0 % use these to extend/shorten spike count window
% motionEnd = 250

a=sum(spikeMatrix(posContrast(posSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
b=sum(spikeMatrix(posContrast(posRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

PCS = mean(a)/.25
PCSerror = sem(a)

PCR = mean(b)/.25
PCRerror = sem(b)

c=sum(spikeMatrix(negContrast(negSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
d=sum(spikeMatrix(negContrast(negRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

NCS = mean(c)/.25
NCSerror=sem(c)

NCR = mean(d)/.25
NCRerror = sem(d)

e=sum(spikeMatrix(posContrast(centerSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
f=sum(spikeMatrix(posContrast(centerRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
g=sum(spikeMatrix(posContrast(center180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

CPCS= mean(e)/.25
CPCSerror = sem(e)
CPCR = mean(f)/.25
CPCRerror = sem(f)
CPC180 = mean(g)/.25
CPC180error = sem(g)

h=sum(spikeMatrix(posContrast(centerRSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
i=sum(spikeMatrix(posContrast(centerRRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
j=sum(spikeMatrix(posContrast(centerR180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

rCPCS= mean(h)/.25
rCPCSerror = sem(h)
rCPCR = mean(i)/.25
rCPCRerror = sem(i)
rCPC180 = mean(j)/.25
rCPC180error = sem(j)

k=sum(spikeMatrix(negContrast(centerSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
l=sum(spikeMatrix(negContrast(centerRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
m=sum(spikeMatrix(negContrast(center180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

CNCS= mean(k)/.25
CNCSerror = sem(k)
CNCR = mean(l)/.25
CNCRerror = sem(l)
CNC180 = mean(m)/.25
CNC180error = sem(m)

n=sum(spikeMatrix(negContrast(centerRSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)
o=sum(spikeMatrix(negContrast(centerRRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
p=sum(spikeMatrix(negContrast(centerR180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)

rCNCS= mean(n)/.2
rCNCSerror = sem(n)
rCNCR = mean(o)/.25
rCNCRerror = sem(o)
rCNC180 = mean(p)/.25
rCNC180error = sem(p)

meanCollection = [PCS PCR NCS NCR CPCS CPCR CPC180 rCPCS rCPCR rCPC180 CNCS CNCR CNC180 rCNCS rCNCR rCNC180];
meanKey = ['PCS ', 'PCR ', 'NCS ', 'NCR ', 'CPCS ', 'CPCR ', 'CPC180 ', 'rCPCS ', 'rCPCR ', 'rCPC180 ', 'CNCS ', 'CNCR ', 'CNC180 ', 'rCNCS ', 'rCNCR ', 'rCNC180 '];
% ((preTime*10)+(motionTime*10)):((preTime*10)+(motionTime*10)+(2500))


subplot(2,3,1)
bar(1:2,[PCS,PCR])
hold on
er = errorbar(1:2,[PCS PCR],[PCSerror PCRerror]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
% errorbar(1:2,PCR,PCRerror)

subplot(2,3,2)
bar(1:2,[NCS,NCR])
hold on
er = errorbar(1:2,[NCS NCR],[NCSerror NCRerror]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off

subplot(2,3,3)
bar(1:3,[CPCS,CPCR,CPC180])
hold on
er = errorbar(1:3,[CPCS CPCR CPC180],[CPCSerror CPCRerror CPC180error]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off

subplot(2,3,4)
bar(1:3,[rCPCS,rCPCR,rCPC180])
hold on
er = errorbar(1:3,[rCPCS rCPCR rCPC180],[rCPCSerror rCPCRerror rCPC180error]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off

subplot(2,3,5)
bar(1:3,[CNCS,CNCR,CNC180])
hold on
er = errorbar(1:3,[CNCS CNCR CNC180],[CNCSerror CNCRerror CNC180error]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off

subplot(2,3,6)
bar(1:3,[rCNCS,rCNCR,rCNC180])
hold on
er = errorbar(1:3,[rCNCS rCNCR rCNC180],[rCNCSerror rCNCRerror rCNC180error]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off

% figure(11)
% 
% title('peaks')
% subplot(2,3,1)
% bar(1:2,[PCSpeak,PCRpeak])
% hold on
% 
% subplot(2,3,2)
% bar(1:2,[NCSpeak,NCRpeak])
% 
% subplot(2,3,3)
% bar(1:3,[CPCSpeak,CPCRpeak,CPC180peak])
% 
% subplot(2,3,4)
% bar(1:3,[rCPCSpeak,rCPCRpeak,rCPC180peak])
% 
% subplot(2,3,5)
% bar(1:3,[CNCSpeak,CNCRpeak,CNC180peak])
% 
% subplot(2,3,6)
% bar(1:3,[rCNCSpeak,rCNCRpeak,rCNC180peak])
% hold off

figure(12)

% subplot(2,3,1)
% plot(mean(psthMatrix(seq,:)),'Color','r')
% hold on
% plot(mean(psthMatrix(rand,:)),'Color','b')
% title('Mean PSTH - All epochs combined')
% legend('seq','rand','location','north')
% title('only separated surround')

subplot(2,3,1)
plot(mean(psthMatrix(posContrast(posSeq),:)),'Color','r')
hold on
plot(mean(psthMatrix(posContrast(posRand),:)),'Color','b')
legend('seq','rand','location','north')
title('positive contrast center, separated surround')

subplot(2,3,2)
plot(mean(psthMatrix(negContrast(negSeq),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(negRand),:)),'Color','b')
legend('seq','rand','location','north')
title('negative contrast center, separated surround')

subplot(2,3,3)
plot(mean(psthMatrix(posContrast(centerSeq),:)),'Color','r','LineWidth',1)
hold on
plot(mean(psthMatrix(posContrast(centerRand),:)),'Color','b','LineWidth',1)
plot(mean(psthMatrix(posContrast(center180),:)),'Color','k','LineWidth',1)
legend('seq','rand','180','location','north')
title('all positive center, sequential surround')

subplot(2,3,4)
plot(mean(psthMatrix(posContrast(centerRSeq),:)),'Color','r','LineWidth',1)
hold on
plot(mean(psthMatrix(posContrast(centerRRand),:)),'Color','b','LineWidth',1)
plot(mean(psthMatrix(posContrast(centerR180),:)),'Color','k','LineWidth',1)
legend('seq','rand','180','location','north')
title('all positive center, random surround')


subplot(2,3,5)
plot(mean(psthMatrix(negContrast(centerSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(center180Neg),:)),'Color','k')
legend('seq','rand','180','location','north')
title('all negative center, sequential surround')

subplot(2,3,6)
plot(mean(psthMatrix(negContrast(centerRSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(centerR180Neg),:)),'Color','k')
legend('seq','rand','180','location','north')
title('all negative center, random surround')

%here for each psth?

figure(13)
for z = 1:length(centerRand)
subplot(1,2,1)
plot(psthMatrix(posContrast(centerSeq(z)),:),'Color','r','LineWidth',.1)
hold on
plot(psthMatrix(posContrast(centerRand(z)),:),'Color','b','LineWidth',.1)
plot(psthMatrix(posContrast(center180(z)),:),'Color','k','LineWidth',.1)
legend('seq','rand','180','location','north')
title('all positive center, sequential surround')


subplot(1,2,2)
plot(psthMatrix(posContrast(centerRSeq(z)),:),'Color','r','LineWidth',.1)
hold on
plot(psthMatrix(posContrast(centerRRand(z)),:),'Color','b','LineWidth',.1)
plot(psthMatrix(posContrast(centerR180(z)),:),'Color','k','LineWidth',.1)
legend('seq','rand','180','location','north')
title('all positive center, random surround')
end






if saveGraph == 1

figure(17)
clf
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
plot(mean(psthMatrix(seq,:)),'Color','r')
hold on
plot(mean(psthMatrix(rand,:)),'Color','b')
plot(mean(psthMatrix(static,:)),'Color','k')
title('Mean PSTH - All epochs combined')
legend('seq','rand','static')

export_fig 'MCSFig.pdf'
    
figure(18)
clf
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
plot(mean(psthMatrix(posContrast(posSeq),:)),'Color','r')
hold on
plot(mean(psthMatrix(posContrast(posRand),:)),'Color','b')
plot(mean(psthMatrix(posContrast(posStatic),:)),'Color','k')
legend('seq','rand','static')
title('positive contrast center - mean PSTH')
export_fig 'MCSFig.pdf' -append


figure(19)
clf
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
plot(mean(psthMatrix(negContrast(negSeq),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(negRand),:)),'Color','b')
plot(mean(psthMatrix(negContrast(negStatic),:)),'Color','k')
legend('seq','rand','static')
title('negative contrast center - mean PSTH')
export_fig 'MCSFig.pdf' -append



figure(20)
clf
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
plot(mean(psthMatrix(posContrast(centerSeq),:)),'Color','r')
hold on
plot(mean(psthMatrix(posContrast(centerRand),:)),'Color','b')
plot(mean(psthMatrix(posContrast(center180),:)),'Color','k')
legend('seq','rand','180')
title('center conditions w/ sequential surround, positive bars - mean PSTH')
export_fig 'MCSFig.pdf' -append




figure(21)
clf
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
plot(mean(psthMatrix(posContrast(centerRSeq),:)),'Color','r')
hold on
plot(mean(psthMatrix(posContrast(centerRRand),:)),'Color','b')
plot(mean(psthMatrix(posContrast(centerR180),:)),'Color','k')
legend('seq','rand','180')
title('center conditions w/ random surround, positive bars - mean PSTH')
export_fig 'MCSFig.pdf' -append



figure(22)
clf
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
plot(mean(psthMatrix(negContrast(centerSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(center180Neg),:)),'Color','k')
legend('seq','rand','180')
title('center conditions w/ sequential surround, negative bars - mean PSTH')
export_fig 'MCSFig.pdf' -append


figure(23)
clf
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
plot(mean(psthMatrix(negContrast(centerRSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(centerR180Neg),:)),'Color','k')
legend('seq','rand','180')
title('center conditions w/ sequential surround, negative bars - mean PSTH')
export_fig 'MCSFig.pdf' -append

end
end

function simpleSpots(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

sampleRate = 10000;
stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;


analysisSplitters = string(splitCell(1,:));
spotIntensityPlace = find(strcmp(analysisSplitters,"spotIntensity")==1);

    for a = 1:size(indexHolder,2)
        spotSizes = splitCell{2,size(splitCell,2)};
        spotSizes = spotSizes((indexHolder{2,a}));
         [spotSizes I] = sort(spotSizes);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);  %sort indices like spotSizes (which are low to high)
         for b = 1:length(unique(spotSizes))
            uniqueSpotSizes = unique(spotSizes); %because I don't know how to index the unique function
            finalInd = find(spotSizes==uniqueSpotSizes(b));    
            onResp(b) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
            offResp(b) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));                
            
            psthData(b,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(b,:)+(100*(b-1)))
            title(indexHolder{1,a})
         end
        
         
         
         figure(1) 
        
         subplot(size(indexHolder,2),1,a)
         if indexHolder{1,a}(spotIntensityPlace) == 0
        plot(uniqueSpotSizes,onResp,'Color','r','LineWidth',2); hold on
        plot(uniqueSpotSizes,offResp,'Color','b','LineWidth',2);
        legend('Off Response','On Response')
        xlabel('spot sizes')
        ylabel('spike count')
        title('Dark Spot')
         else
        plot(uniqueSpotSizes,onResp,'Color','b','LineWidth',2); hold on
        plot(uniqueSpotSizes,offResp,'Color','r','LineWidth',2);
        legend('On Response','Off Response')
        xlabel('spot sizes')
        ylabel('spike count')
        
        title("Light Spot --" + indexHolder{1,a}(2)) %will need to automate later
         end
      
    end
         

end

function [MNOFFParasol,MNONParasol]=CRF(spikeMatrix,psthMatrix,radius,timings,saveStuff,MNOFFParasol,MNONParasol)


timings = timings * 10;
uniqueR = unique(radius);
uniqueR = sort(uniqueR);
xaxis = [];
yaxis = [];
for p = 1:length(uniqueR)
    
    cSpot = find(radius == uniqueR(p));
    xaxis = [xaxis mean(radius(cSpot))];
    yaxis = [yaxis mean(sum(spikeMatrix(cSpot',timings(1):timings(1)+timings(2)),2))];

end

xaxis
yaxis

if saveStuff.flag == 1 && strcmp(saveStuff.cellType,'OFF Parasol')
    disp('accessed')
    MNOFFParasol.Exp(saveStuff.expNum).Cell(saveStuff.OFFParasolC).Sensitivity = [xaxis;yaxis];
elseif saveStuff.flag == 1 && strcmp(saveStuff.cellType,'ON Parasol')
    MNONParasol.Exp(saveStuff.expNum).Cell(saveStuff.ONParasolC).Sensitivity = [xaxis;yaxis];
end
plot(xaxis,yaxis)
axis([0 max(xaxis) 0 max(yaxis)+2])

end

function mTF(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)
sampleRate = 10000;

stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;


analysisSplitters = string(splitCell(1,:));
tFreq= find(strcmp(analysisSplitters,"temporalFrequency")==1);
figure(10); clf; hold on


    for a = 1:size(indexHolder,2)
        radii = splitCell{2,size(splitCell,2)};
        radii = radii((indexHolder{2,a}));
         [radii I] = sort(radii);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I)
         
         

        for u = 1:length(unique(radii))
            uniqueRadii = unique(radii); %because I don't know how to index the unique function
            finalInd = find(radii==uniqueRadii(u));
%             angleIndex = find(orientation == uniqueAngle(u));
%             masterIndex = intersect(angleIndex,apertureIndex(a,:));
            data = mean(spikeMatrix(sortedIndex(finalInd),:),1);
                    binnedData = BinSpikeRate(data(stimStart:stimEnd), 100, sampleRate);
                [F, phase] = frequencyModulation(binnedData, ...
                100, indexHolder{1,a}(1,params.tfreq(1)), 'avg', 1:2, []);
                avgF1(u)= F(1);
                avgF2(u)= F(2);    
                
                phaser(u,:) = phase;
                spikeCount(u) = sum(data(stimStart:stimEnd),2);
                
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
        end
        
%         uniqueRadii

    figure(a)
    subplot(3,1,1);
    plot(uniqueRadii, avgF1(1:length(uniqueRadii))*1e-4);
%     titleString = strcat({'F1 Resp'}, num2str(indexHolder{1,1}));
    title('F1 Resp')
    xlabel(analysisSplitters(1:3))
    hold on
    subplot(3,1,2);
    plot(uniqueRadii, avgF2(1:length(uniqueRadii))*1e-4);
%     titleString = strcat({'F2 Resp'}, num2str(indexHolder{1,1}));
    title('F2 Resp')
    xlabel(num2str(indexHolder{1,a}'))
    
    subplot(3,1,3);
    plot(uniqueRadii, phaser(1:length(uniqueRadii))*1e-4);
%     titleString = strcat({'Spike Count'}, num2str(indexHolder{1,1}));
    title('p h a s e')
    xlabel(analysisSplitters(1:3))
    
    figure(10)

   
    
    end

end

function old(spikeMatrix,psthMatrix,radius,timings,tFrequency,sampleRate,whichStim)



sender = struct();
uniqueR = unique(radius);
uniqueR = sort(uniqueR);
% stimIndex = find(whichStim == 1);
for s = 1:length(uniqueR)
spotIndex = find(radius == uniqueR(s));
spotIndex = spotIndex';
if isempty(find(whichStim==1))
stimIndex = find(whichStim(spotIndex) == 2);
else
stimIndex = find(whichStim(spotIndex) == 1);   
end
stimIndex
sender(s).data = mean(spikeMatrix(spotIndex(stimIndex),:),1);
% size(mean(spikeMatrix(spotIndex,:)));
sender(s).params.spatialFrequency = uniqueR(s);
sender(s).params.temporalFrequency = tFrequency;
sender(s).params.preTime = timings(1)*1e-3;
sender(s).params.stimStart = (timings(1)*1e-3)*sampleRate+1;
sender(s).params.stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
sender(s).params.sampleRate = sampleRate;
end
result = sMTFAnalysis(sender,1e4,'spikes','avg');
figure(1)
    hold on
    plot(result.uniqueSF,result.avgF1);
    plot(result.uniqueSF,result.avgF2);
    title('F1 Response to Spot')
    xlabel('Spot Radius (micron)')
    ylabel('F1 Amplitude (Hz)')
end

function centeringBars(XspikeMatrix,YspikeMatrix,positionX,positionY,timings,tFrequency,sampleRate)

sender = struct();
uniqueP = unique(positionX);
uniqueP = sort(uniqueP);
for s = 1:length(uniqueP)
barIndex = find(positionX == uniqueP(s));
barIndex = barIndex';
sender(s).data = mean(XspikeMatrix(barIndex,:),1);
sender(s).params.spatialFrequency = uniqueP(s);
sender(s).params.temporalFrequency = tFrequency;
sender(s).params.preTime = timings(1)*1e-3;
sender(s).params.stimStart = (timings(1)*1e-3)*sampleRate+1;
sender(s).params.stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
sender(s).params.sampleRate = sampleRate;
end
result = sMTFAnalysis(sender,1e4,'spikes','avg');
figure(1)
subplot(1,2,1)
    plot(result.uniqueSF,result.avgF1);
    hold on
    plot(result.uniqueSF,result.avgF2);
    title('F1 Response to Spot')
    xlabel('Spot Radius (micron)')
    ylabel('F1 Amplitude (Hz)')

uniqueP = unique(positionY);
uniqueP = sort(uniqueP);
    
for s = 1:length(uniqueP)
    barIndex = find(positionY == uniqueP(s));
    barIndex = barIndex';
    sender(s).data = mean(YspikeMatrix(barIndex,:),1);

    sender(s).params.spatialFrequency = uniqueP(s);
    sender(s).params.temporalFrequency = tFrequency;
    sender(s).params.preTime = timings(1)*1e-3;
    sender(s).params.stimStart = (timings(1)*1e-3)*sampleRate+1;
    sender(s).params.stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
    sender(s).params.sampleRate = sampleRate;
end
result = sMTFAnalysis(sender,1e4,'spikes','avg');
subplot(1,2,2)
    plot(result.uniqueSF,result.avgF1);
    hold on
    plot(result.uniqueSF,result.avgF2);
    title('F1 Response to Spot')
    xlabel('Spot Radius (micron)')
    ylabel('F1 Amplitude (Hz)')

end

function orientedStim(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

% for z = 1:length(splitCell)-1    
%     uniqueSplit = unique(splitCell{2,z}); 
%     for c = 1:length(uniqueSplit)
%     splitIndex = find(splitCell{2,z}==uniqueSplit(c));
%     
%     end
% end
sampleRate = 10000;
stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;


analysisSplitters = string(splitCell(1,:));
tFreq= find(strcmp(analysisSplitters,"temporalFrequency")==1);
figure(10); clf; hold on
if strcmp(params.stimName, 'Grating')

    for a = 1:size(indexHolder,2)
        orientations = splitCell{2,size(splitCell,2)};
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
         
         

        for u = 1:length(unique(orientations))
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
%             angleIndex = find(orientation == uniqueAngle(u));
%             masterIndex = intersect(angleIndex,apertureIndex(a,:));
            data = mean(spikeMatrix(sortedIndex(finalInd),:),1);
                    binnedData = BinSpikeRate(data(stimStart:stimEnd), 100, sampleRate);
                [F, phase] = frequencyModulation(binnedData, ...
                100, indexHolder{1,a}(1,tFreq), 'avg', 1:2, []);
                avgF1(u)= F(1);
                avgF2(u)= F(2);     
                spikeCount(u) = sum(data(stimStart:stimEnd),2);
                
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(15); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
        end

    figure(a)
    subplot(1,3,1);
    polar(uniqueOrientations * 2 * pi/360, avgF1(1:length(uniqueOrientations))*1e-4);
%     titleString = strcat({'F1 Resp'}, num2str(indexHolder{1,1}));
    title('F1 Resp')
    xlabel(analysisSplitters(1:3))
    hold on
    subplot(1,3,2);
    polar(uniqueOrientations * 2 * pi/360, avgF2(1:length(uniqueOrientations))*1e-4);
%     titleString = strcat({'F2 Resp'}, num2str(indexHolder{1,1}));
    title('F2 Resp')
    xlabel(num2str(indexHolder{1,a}'))
    
    subplot(1,3,3);
    polar(uniqueOrientations * 2 * pi/360, spikeCount(1:length(uniqueOrientations)));
%     titleString = strcat({'Spike Count'}, num2str(indexHolder{1,1}));
    title('Spike Count')
    xlabel(analysisSplitters(1:3))
    
    figure(10)

   
    
    end

elseif strcmp(params.stimName,'bars')
    disp('bars')


figure(10);clf;
onResp =[];
offResp=[];
for a = 1:size(indexHolder,2)
        orientations = splitCell{2,size(splitCell,2)};
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
        clear onResp
        clear offResp
    
    for u = 1:length(unique(orientations))
        
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
            onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
            offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));
            
            %plot traces
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(15); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
    end



% size(onResp)
% uniqueOrientations = unique(splitCell{2,size(splitCell,2)});
% uniqueOrientations
% onResp



figure(a)
subplot(1,2,1);
polar(uniqueOrientations * 2 * pi /360, onResp);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
hold on
subplot(1,2,2);
polar(uniqueOrientations * 2 * pi/360, offResp);
title('OFF')
xlabel(num2str(indexHolder{1,a}'));

end
elseif strcmp(params.stimName,'moving bar')
    figure(10);clf;
resp = [];
for a = 1:size(indexHolder,2)
        orientations = splitCell{2,size(splitCell,2)};
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
        clear resp
    
    for u = 1:length(unique(orientations))
        
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
            resp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),:),2));
            
            
            %plot traces
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
            
    end

resp
uniqueOrientations    

% if isequal(size(uniqueOrientations),size(resp))
%     disp('equal')
% holdOrientations = uniqueOrientations;
% elseif ~isequal(size(uniqueOrientations),size(resp))
% uniqueOrientations = holdOrientations;
% end

size(uniqueOrientations)
size(resp)


figure(a)
polar(uniqueOrientations * 2 * pi /360, resp);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
ylabel(num2str(indexHolder{1,a}'));




end



end
end