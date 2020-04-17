%% initialize save structs
MNOFFParasol = struct();
MNONParasol = struct();
MNONSmooth = struct();
MNOFFSmooth = struct();

% MCS -- need something in case no cell

CMSOFFParasol = [];
CMSONParasol = [];
CMSONSmooth = [];
CMSOFFSmooth = [];
%% I'M DUMB
OFFParasolC = 0;
ONParasolC = 0;
ONSmoothC = 0;
OFFSmoothC = 0;

OFFParasolC2 =0;
ONParasolC2 = 0;
ONSmoothC2 = 0;
OFFSmoothC2 = 0;
%%

% SAVE PARAMS HERE:
cellType = 'ON Parasol';
saveFlag = 0;
lastRun = 0;
expNum= 11;
cellName = 'Bc8';

% load FT first, give different name --- #1
cd('/Users/toddappleby/Documents/Data/Clarinet Exports/2020_0318')

load('20200318Bc4_FT.mat')
frameTimings = epochs;
expDate = dir;
% now load data -- #2
load('20200318Bc4.mat')

uniqueProtocols = [];

for z = 1:size(epochs,2)
   list(z) = string(epochs(z).meta.displayName);
%    allEpochData(z,1:length(epochs(z).epoch)) = epochs(z).epoch;
end
%find 0s, make index, create new string
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
   egLabel = epochs(i).meta.epochGroupLabel;
   
   if strcmp(displayName,'Motion And Noise') && strcmp(egLabel,'Control')
        
        
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
 frameDwell(count) = epochs(i).meta.frameDwell;
% frameDwell(count)=1;
noiseClass(count) = {epochs(i).meta.noiseClass};
barWidth(count) = epochs(i).meta.barWidth;
apRadius(count) = epochs(i).meta.apertureRadius;
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

desiredSTD = 5; %arbitrary, works for most extracellular data

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
frameDwell = frameDwell(DI2);
 barWidth = barWidth(DI2);
 apRadius = apRadius(DI2);
 bgClass = bgClass(DI2);
 seedList = seedList(DI2);
 noiseClass = noiseClass(DI2);
%  binaryIndex = find(contains(noiseClass,'binary'));
%  gaussianIndex = find(contains(noiseClass,'gaussian'));
binaryIndex = contains(noiseClass,'binary');
gaussianIndex = contains(noiseClass,'gaussian');
 bgClassOrig = bgClass;
 %% split by protocol variable...somehow
 %run unique(apRadius) or unique(barWidth)
%  splitter = 300;
% splitIndex = find(apRadius == splitter);

% seedList = seedList(splitIndex);
% spikeMatrix = spikeMatrix(splitIndex,:);
% psthMatrix = psthMatrix(splitIndex,:);
% bgClass = bgClass(splitIndex);
% gaussianIndex = gaussianIndex(splitIndex);
% binaryIndex = binaryIndex(splitIndex);
 
% NOTE: THIS DOES NOT WORK IF GAUSSIAN AND BINARY ARE BEING SORTED BY SAME
% FRAME DWELL (you'll get logical 1s in binaryIndex where it's actually
% gaussian noise)

splitter2 = 1;




splitIndex2 = frameDwell==splitter2;
% bgClass = bgClass(splitIndex2>0);
binaryIndex = splitIndex2==1;     
% frameDwell = frameDwell(splitIndex2);
 
 %% Get stim frames -- #5 (standard, using Mike's frame time functions)
count = 0;
noiseFlag =0; %1 for gaussian
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
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,1,seedList(gaussianIndex),frameDwell(gaussianIndex));
else
    noiseVars = struct();
noiseVars.type = 'binary';
noiseVars.contrast = 0.3333;

timings = [250,10000,250]; % AUTOMATE THIS LATER
frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage(binaryIndex,:), 10000, binRate);
% frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage, 10000, 10000);
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,1,seedList(binaryIndex),frameDwell(binaryIndex));
end

% just use this for gaussian/binary - #6 (keeps 1 seeds for now)
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
 
% correlation stuff -- #7
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


% seqPSTH(countSeq,:) = psthMatrix(h,:);

    elseif bgClass(h) == 2
        countR = countR+1;
lfResp = response1(h,:);
lfStim = frameValues(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterRandom = (linearFilterRandom*(countR-1) + lf)/countR;
sortNLRand(countR) = h;

% randPSTH(countR,:) = psthMatrix(h,:);

    elseif bgClass(h) == 3
        countStatic = countStatic + 1;
lfResp = response1(h,:);
lfStim = frameValues(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterStatic = (linearFilterStatic*(countStatic-1) + lf)/countStatic;
sortNLStatic(countStatic) = h;


% staticPSTH(countStatic,:) = psthMatrix(h,:);

    end
        
end

% figure(3)
% subplot(3,1,1)
% plot(mean(seqPSTH)*10000,'r')
% subplot(3,1,2)
% plot(mean(randPSTH)*10000,'b')
% subplot(3,1,3)
% plot(mean(staticPSTH)*10000,'k')
% sgtitle('Mean PSTH over all Gaussian Trials')
% xlabel('Time (ms)')
% ylabel('Spike Rate (Hz)')
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')

% figure(2)
% plot(mean(seqPSTH)*10000,'r')
% hold on
% plot(mean(randPSTH)*10000)

% plot(mean(staticPSTH)*10000,'k')
% graphX = linspace(0,10,10500);
% cheese = table((mean(seqPSTH)*10000)',(mean(randPSTH)*10000)',(mean(staticPSTH)*10000)',graphX');
% writetable(cheese,'thing.xlsx','sheet',1);



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
     

% plotlf - raw filters #8

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

% filterTable = table(linearFilterSeq(1:plotLngth)',linearFilterRandom(1:plotLngth)',linearFilterStatic(1:plotLngth)',((1:plotLngth)/binRate)');
% writetable(filterTable,'FilterTable.xlsx','sheet',1);

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
 
 % start NL 
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
 %max(nlX)
 plot(nlX,5000*outputNonlinearity(nlParams,nlX),'Color','r','LineWidth',2)
 hold on
 
 sequentialOut = 5000*outputNonlinearity(nlParams,nlX);
 
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
 
% random NL 

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
 %max(nlX)
 plot(nlX,5000*outputNonlinearity(nlParams,nlX),'Color','b','LineWidth',2)
 hold on
 randomOut = 5000*outputNonlinearity(nlParams,nlX);
 
 if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
   MNOFFParasol.Exp(expNum).Cell(OFFParasolC).NL.Rand = outputNonlinearity(nlParams,nlX);
   elseif saveFlag ==1 && strcmp(cellType,'ON Parasol')
   MNONParasol.Exp(expNum).Cell(ONParasolC).NL.Rand = outputNonlinearity(nlParams,nlX);
   elseif saveFlag ==1 && strcmp(cellType,'ON Smooth')
   
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).NL.Rand = outputNonlinearity(nlParams,nlX);
   elseif saveFlag ==1 && strcmp(cellType,'OFF Smooth')
    
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).NL.Rand = outputNonlinearity(nlParams,nlX);
 end
% stationary NL 

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
 %max(nlX)
 figure(9)
 plot(nlX,5000*outputNonlinearity(nlParams,nlX),'Color','k','LineWidth',2)
  legend('sequential','random','static','location','northwest')
  title('Non Linearity for each surround condition')
  xlabel('input','FontSize',18)
  ylabel('spikes/s','FontSize',18)
  set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
export_fig 'noiseNL2.pdf' -append

if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
    MNOFFParasol.Exp(expNum).Cell(OFFParasolC).NL.Static = outputNonlinearity(nlParams,nlX);
elseif saveFlag == 1 && strcmp(cellType,'ON Parasol')
    MNONParasol.Exp(expNum).Cell(ONParasolC).NL.Static = outputNonlinearity(nlParams,nlX);
elseif saveFlag ==1 && strcmp(cellType,'ON Smooth')
    MNONSmooth.Exp(expNum).Cell(ONSmoothC).NL.Static = outputNonlinearity(nlParams,nlX);
elseif saveFlag ==1 && strcmp(cellType,'OFF Smooth')
    MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).NL.Static = outputNonlinearity(nlParams,nlX);
end

xAxisNL = nlX/max(nlX);
staticOut = 5000*outputNonlinearity(nlParams,nlX);

outTable = table(xAxisNL', staticOut', randomOut',sequentialOut');
writetable(outTable,'outTable2.xlsx','sheet',1);


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
save('MotionandNoise.mat','-append','MNONSmooth')
save('MotionandNoise.mat','-append','MNOFFSmooth')

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
    OFFParasolC2 = OFFParasolC2 + 1;
CMSOFFParasol(OFFParasolC2,:) = meanCollection;
    elseif saveFlag == 1 && strcmp(cellType,'ON Parasol')
            ONParasolC2 = ONParasolC2 + 1;
            CMSONParasol(ONParasolC2,:) = meanCollection;
        elseif saveFlag == 1 && strcmp(cellType,'ON Smooth')
                ONSmoothC2 = ONSmoothC2 + 1;
                CMSONSmooth(ONSmoothC2,:) = meanCollection;
            elseif saveFlag == 1 && strcmp(cellType,'OFF Smooth')
                    OFFSmoothC2 = OFFSmoothC2 + 1;
                    CMSOFFSmooth(OFFSmoothC2,:) = meanCollection;
end


if saveFlag == 1 && lastRun == 1
    firstSlot=3;
    if strcmp(expDate(3).name(1:3),'.DS')
        firstSlot = 4;
    end
testSave = strcat(expDate(firstSlot).name(1:8),'_CMS');
    cd('/Users/toddappleby/Documents/Data/Clarinet Exports/SavedData')
    
    
    save(testSave,'CMSOFFParasol')
    
    save(testSave,'-append','CMSONParasol')
    save(testSave,'-append','CMSONSmooth')
    save(testSave,'-append','CMSOFFSmooth')
end


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
desiredSTD = 5;
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

%% functions!!
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

function [meanCollection, meanKey] = motionCS(epochs,bgClass,centerClass,contrastState,epochStorage,protocolExample,saveGraph)

sampleRate = 10000;
binRate = 1000;
preTime = epochs(protocolExample).meta.preTime;
stimTime = epochs(protocolExample).meta.stimTime;
stimOrig = stimTime;
motionTime = epochs(protocolExample).meta.motionTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = 5.5;
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


 
posContrast = find(contrastState==abs(unique(contrastState(1))));
 negContrast = find(contrastState==-abs(unique(contrastState(1))));

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


% returnStruct = struct();
% returnStruct.surround.PSeq = mean(psthMatrix(posContrast(posSeq),:));
% returnStruct.surround.PRand = mean(psthMatrix(posContrast(posRand),:));
% returnStruct.surround.NSeq = mean(psthMatrix(negContrast(negSeq),:));
% returnStruct.surround.NRand = mean(psthMatrix(negContrast(negRand),:));
% returnStruct.center.PCS = mean(psthMatrix(posContrast(centerSeq),:));
% returnStruct.center.PCR = mean(psthMatrix(posContrast(centerRand),:));
% returnStruct.center.PCSt = mean(psthMatrix(posContrast(center180),:));
% returnStruct.center.rPCS = mean(psthMatrix(posContrast(centerRSeq),:));
% returnStruct.center.rPCR = mean(psthMatrix(posContrast(centerRRand),:));
% returnStruct.center.rPCSt = mean(psthMatrix(posContrast(centerR180),:));
% returnStruct.center.NCS = mean(psthMatrix(negContrast(centerSeqNeg),:));
% returnStruct.center.NCR = mean(psthMatrix(negContrast(centerRandNeg),:));
% returnStruct.center.NCSt = mean(psthMatrix(negContrast(center180Neg),:));
% returnStruct.center.rNCS = mean(psthMatrix(negContrast(centerRSeqNeg),:));
% returnStruct.center.rNCR = mean(psthMatrix(negContrast(centerRRandNeg),:));
% returnStruct.center.rNCSt = mean(psthMatrix(negContrast(centerR180Neg),:));



figure(10)
title('means')
preTime = preTime/10;
motionTime=motionTime/10;

motionStart = 100; % use these to extend/shorten spike count window
motionEnd = 250;
mDuration=100;
% 
% psthMatrix(negContrast(negSeq(1)),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd)))
% psthMatrix(negContrast(negSeq(1)),:)
% 
% %  a=mean(psthMatrix(posContrast(posSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)
% % b=mean(psthMatrix(posContrast(posRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)
% psthMatrix(posContrast(posSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd)))
% a=sum(psthMatrix(posContrast(posSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% b=sum(psthMatrix(posContrast(posRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% 
% PCS = mean(a);
% PCSpeak = max(a);
% PCSerror = sem(a);
% 
% PCR = mean(b);
% PCRpeak = max(b);
% PCRerror = sem(b);
% 
% c=mean(sum(psthMatrix(negContrast(negSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2))/mDuration;
% d=mean(sum(psthMatrix(negContrast(negRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2))/mDuration;
% 
% NCS = mean(c);
% NCSpeak = max(c);
% NCSerror=sem(c);
% 
% NCR = mean(d);
% NCRpeak = max(d);
% NCRerror = sem(d);
% 
% e=sum(psthMatrix(posContrast(centerSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% f=sum(psthMatrix(posContrast(centerRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% g=sum(psthMatrix(posContrast(center180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% 
% CPCS= mean(e);
% CPCSpeak = max(e);
% CPCSerror = sem(e);
% CPCR = mean(f);
% CPCRpeak = max(f);
% CPCRerror = sem(f);
% CPC180 = mean(g);
% CPC180peak = max(g);
% CPC180error = sem(g);
% 
% h=sum(psthMatrix(posContrast(centerRSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% i=sum(psthMatrix(posContrast(centerRRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% j=sum(psthMatrix(posContrast(centerR180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% 
% rCPCS= mean(h);
% rCPCSpeak = max(h);
% rCPCSerror = sem(h);
% rCPCR = mean(i);
% rCPCRpeak = max(i);
% rCPCRerror = sem(i);
% rCPC180 = mean(j);
% rCPC180peak = max(j);
% rCPC180error = sem(j);
% 
% k=sum(psthMatrix(negContrast(centerSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% l=sum(psthMatrix(negContrast(centerRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% m=sum(psthMatrix(negContrast(center180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% 
% CNCS= mean(k);
% CNCSpeak = max(k);
% CNCSerror = sem(k);
% CNCR = mean(l);
% CNCRpeak = max(l);
% CNCRerror = sem(l);
% CNC180 = mean(m);
% CNC180peak = max(m);
% CNC180error = sem(m);
% 
% n=sum(psthMatrix(negContrast(centerRSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% o=sum(psthMatrix(negContrast(centerRRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% p=sum(psthMatrix(negContrast(centerR180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2)/mDuration;
% 
% rCNCS= mean(n);
% rCNCSpeak = max(n);
% rCNCSerror = sem(n);
% rCNCR = mean(o);
% rCNCRpeak = max(o);
% rCNCRerror = sem(o);
% rCNC180 = mean(p);
% rCNC180peak=max(p);
% rCNC180error = sem(p);


% 
% 
% 
size(spikeMatrix);

a=sum(spikeMatrix(posContrast(posSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
b=sum(spikeMatrix(posContrast(posRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
a=a/.1;
b=b/.1;


PCS = mean(a);
PCSerror = sem(a);

PCR = mean(b);
PCRerror = sem(b);

c=sum(spikeMatrix(negContrast(negSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
d=sum(spikeMatrix(negContrast(negRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
c=c/.1;
d=d/.1;



NCS = mean(c);
NCSerror=sem(c);

NCR = mean(d);
NCRerror = sem(d);

e=sum(spikeMatrix(posContrast(centerSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
f=sum(spikeMatrix(posContrast(centerRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
g=sum(spikeMatrix(posContrast(center180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

e=e/.1;
f=f/.1;
g=g/.1;


CPCS= mean(e);
CPCSerror = sem(e);
CPCR = mean(f);
CPCRerror = sem(f);
CPC180 = mean(g);
CPC180error = sem(g);

h=sum(spikeMatrix(posContrast(centerRSeq),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
i=sum(spikeMatrix(posContrast(centerRRand),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
j=sum(spikeMatrix(posContrast(centerR180),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

h=h/.1;
i=i/.1;
j=j/.1;

rCPCS= mean(h);
rCPCSerror = sem(h);
rCPCR = mean(i);
rCPCRerror = sem(i);
rCPC180 = mean(j);
rCPC180error = sem(j);

k=sum(spikeMatrix(negContrast(centerSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
l=sum(spikeMatrix(negContrast(centerRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
m=sum(spikeMatrix(negContrast(center180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

k=k/.1;
l=l/.1;
m=m/.1;

CNCS= mean(k);
CNCSerror = sem(k);
CNCR = mean(l);
CNCRerror = sem(l);
CNC180 = mean(m);
CNC180error = sem(m);

n=sum(spikeMatrix(negContrast(centerRSeqNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
o=sum(spikeMatrix(negContrast(centerRRandNeg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);
p=sum(spikeMatrix(negContrast(centerR180Neg),((preTime*10)+(motionTime*10)-motionStart):((preTime*10)+(motionTime*10)+(motionEnd))),2);

n=n/.1;
o=o/.1;
p=p/.1;

rCNCS= mean(n);
rCNCSerror = sem(n);
rCNCR = mean(o);
rCNCRerror = sem(o);
rCNC180 = mean(p);
rCNC180error = sem(p);

meanCollection = [PCS PCR NCS NCR CPCS CPCR CPC180 rCPCS rCPCR rCPC180 CNCS CNCR CNC180 rCNCS rCNCR rCNC180];
meanKey = ['PCS ', 'PCR ', 'NCS ', 'NCR ', 'CPCS ', 'CPCR ', 'CPC180 ', 'rCPCS ', 'rCPCR ', 'rCPC180 ', 'CNCS ', 'CNCR ', 'CNC180 ', 'rCNCS ', 'rCNCR ', 'rCNC180 '];



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
legend('seq','rand','location','northwest')
title('positive contrast center, separated surround')

subplot(2,3,2)
plot(mean(psthMatrix(negContrast(negSeq),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(negRand),:)),'Color','b')
legend('seq','rand','location','northwest')
title('negative contrast center, separated surround')

subplot(2,3,3)
plot(mean(psthMatrix(posContrast(centerSeq),:)),'Color','r','LineWidth',1)
hold on
plot(mean(psthMatrix(posContrast(centerRand),:)),'Color','b','LineWidth',1)
plot(mean(psthMatrix(posContrast(center180),:)),'Color','k','LineWidth',1)
legend('seq','rand','180','location','northwest')
title('all positive center, sequential surround')

subplot(2,3,4)
plot(mean(psthMatrix(posContrast(centerRSeq),:)),'Color','r','LineWidth',1)
hold on
plot(mean(psthMatrix(posContrast(centerRRand),:)),'Color','b','LineWidth',1)
plot(mean(psthMatrix(posContrast(centerR180),:)),'Color','k','LineWidth',1)
legend('seq','rand','180','location','northwest')
title('all positive center, random surround')


subplot(2,3,5)
plot(mean(psthMatrix(negContrast(centerSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(center180Neg),:)),'Color','k')
legend('seq','rand','180','location','northwest')
title('all negative center, sequential surround')

subplot(2,3,6)
plot(mean(psthMatrix(negContrast(centerRSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(centerR180Neg),:)),'Color','k')
legend('seq','rand','180','location','northwest')
title('all negative center, random surround')

figure(50)
subplot(1,2,1)
plot(mean(psthMatrix(negContrast(centerSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(center180Neg),:)),'Color','k')
ylabel('Spike/sec','FontSize',14)
xlabel('Time (ms)','Fontsize',14)
legend('Motion','Random','Motion-Reversed','location','northwest')
title('3 Center Conditions w/ Surround Motion','FontSize',12)


subplot(1,2,2)
plot(mean(psthMatrix(negContrast(centerRSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(centerR180Neg),:)),'Color','k')
ylabel('Spike/sec','FontSize',14)
xlabel('Time (ms)','Fontsize',14)
legend('Motion','Random','Motion-Reversed','location','northwest')
title('3 Center Conditions w/ Non-motion Surround','FontSize',12)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
export_fig 'noiseNL2.pdf' -append

figure(51)
subplot(1,2,1)
plot(mean(psthMatrix(posContrast(centerSeq),:)),'Color','r','LineWidth',1)
hold on
plot(mean(psthMatrix(posContrast(centerRand),:)),'Color','b','LineWidth',1)
plot(mean(psthMatrix(posContrast(center180),:)),'Color','k','LineWidth',1)
legend('seq','rand','180','location','northwest')
title('all positive center, sequential surround')

subplot(1,2,2)
plot(mean(psthMatrix(posContrast(centerRSeq),:)),'Color','r','LineWidth',1)
hold on
plot(mean(psthMatrix(posContrast(centerRRand),:)),'Color','b','LineWidth',1)
plot(mean(psthMatrix(posContrast(centerR180),:)),'Color','k','LineWidth',1)
legend('seq','rand','180','location','northwest')
title('all positive center, random surround')

%here for each psth?

% figure(13)
% for z = 1:length(centerRand)
% subplot(1,2,1)
% plot(psthMatrix(posContrast(centerSeq(z)),:),'Color','r','LineWidth',.1)
% hold on
% plot(psthMatrix(posContrast(centerRand(z)),:),'Color','b','LineWidth',.1)
% plot(psthMatrix(posContrast(center180(z)),:),'Color','k','LineWidth',.1)
% legend('seq','rand','180','location','north')
% title('all positive center, sequential surround')
% 
% 
% subplot(1,2,2)
% plot(psthMatrix(posContrast(centerRSeq(z)),:),'Color','r','LineWidth',.1)
% hold on
% plot(psthMatrix(posContrast(centerRRand(z)),:),'Color','b','LineWidth',.1)
% plot(psthMatrix(posContrast(centerR180(z)),:),'Color','k','LineWidth',.1)
% legend('seq','rand','180','location','north')
% title('all positive center, random surround')
% end


end