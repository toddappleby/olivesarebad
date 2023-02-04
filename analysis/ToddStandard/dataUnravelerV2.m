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
cellName = 'Ac1';

% % load FT first, give different name --- #1
% cd('/Users/toddappleby/Documents/Data/Clarinet Exports/2020_0630')
% %smooth
% cd('/Users/toddappleby/Documents/Data/Clarinet Exports/2021_0128') %parasol
% cd('/Users/toddappleby/Documents/Data/Clarinet Exports/2020_0713')%BT?
% cd('/Users/toddappleby/Documents/Data/Clarinet Exports/2021_0225')
% cd('E:\Data Analysis_2020\2020_0519')

% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2019_1010')

% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2020_0504')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2021_0128')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2021_0507')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2021_0812')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2020_0206')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2021_0427')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2020_0713')
cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2020_0630')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2021_0910')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2021_0907')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2021_1102')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2020_0930')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2020_0713')

% load('20200713Bc2_FT.mat')
% load('20200504Bc1_FT.mat')
% load('20210427Ac1_FT.mat')
% load('20201118Bc3_FT.mat')
% load('20210507Bc4_FT.mat')
% load('20200923Bc2_FT.mat')
% load('20210507Bc4_FT.mat')
% load('20210122Ac7_FT.mat')
load('20200630Ac1_FT.mat')
% load('20210910Fc1_FT.mat')
% load('20210907Ac2_FT.mat')
% load('20210128Bc2_FT.mat')
% load('20200713Bc3_FT.mat')
% load('20200930Bc2_FT.mat')
% load('20211102Ac2_FT.mat')
% load('20191126Bc2_FT.mat')
frameTimings = epochs;
expDate = dir;
% now load data -- #2+
% load('20200713Bc2.mat')
% load('20200504Bc1.mat')
% % load('20201118Bc3.mat')
% load('20210910Fc1.mat')
% load('20200923Bc2.mat')
% load('20210907Ac2.mat')
% load('20200713Bc3.mat')
load('20200630Ac1.mat')
% load('20210128Bc2.mat')
% load('20200825Ac3.mat')
% % load('20210122Ac7.mat')
% load('20200930Bc2.mat')
% load('20211102Ac2.mat')
% load('20191126Bc2.mat')
% load('20210128Bc2.mat')
% load('20210427Ac1.mat')
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
% clear epochNumFrames
for f = 1:length(frameTimings)

     displayName = frameTimings(f).meta.displayName;
    
     
   
   if strcmp(displayName,'Motion And Noise')
       
    count = count+1;
    singleMonitorRun = frameTimings(f).epoch;
 monitorStorage(count,:) = singleMonitorRun;
 epochStartTime(count,:) = frameTimings(f).meta.epochStartTime;
%  epochNumFrames(count) = frameTimings(f).meta.epochNum;
   end
end

[dates,DI] = sort(epochStartTime);  
monitorStorage = monitorStorage(DI,:);
% epochNumFrames = epochNumFrames(DI);
%%  spike detection and epoch organization -- #4
dataType=1;
excorinh = 'Whole cell_exc';
 count = 0;
 desiredSTD =5; 
clear binnedStorage
clear epochStorage
clear epochStartTimeD

for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
  if dataType == 1
%       dTypeLabel = 'Control';
      dTypeLabel = 'motion and noise';
  else
      dTypeLabel = excorinh;
  end
         if strcmp(displayName,'Motion And Noise') && strcmp(egLabel,dTypeLabel)
 
        
        
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
 epochStorage(count,:) = singleEpoch;
%  epochNumResponses(count) = epochs(i).meta.epochNum;
 
seedList(count) = epochs(i).meta.seed;

% epochStartTimeD(count) = epochs(i).meta.epochTime;
epochStartTimeD(count,:) = epochs(i).meta.epochStartTime;
if isfield(epochs(i).meta,'frameDwell')
frameDwell(count) = epochs(i).meta.frameDwell;
else
    frameDwell(count)=1;
end
noiseClass(count) = {epochs(i).meta.noiseClass};
barWidth(count) = epochs(i).meta.barWidth;
apRadius(count) = epochs(i).meta.apertureRadius;
barOrientation(count) = epochs(i).meta.barOrientation;
barFrameDwell(count) = epochs(i).meta.barFrameDwell;
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



% binRate = 1000;
% response(o,:) = binSpikeCount(epoch(o,:)/sampleRate, binRate, sampleRate);
% response(o,:) = psth(response(o,:)*binRate,6+2/3,binRate,1);

%spike detection
if dataType==1
for k = 1:size(binnedStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(binnedStorage(k,:), stimTime, ...
              desiredSTD);
          
          if isempty(spikes)
              disp('deleted epoch')
          else
              
        spikeMatrixUnbinned(k,:) = spikes;
%         psthMatrix2(k,:) = psth(spikeMatrix(k,:)*10000,6+2/3,10000,1);
        spikes = binSpikeCount(spikes, binRate, sampleRate);
        spikes = spikes';
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:)*binRate,20,binRate,1);
%         
%         spikeMatrix(k,:) = binSpikeCount(spikes/sampleRate, binRate, sampleRate);
%         psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
        
          end
end

else
    clear spikeMatrix
    clear psthMatrix

 for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:250));
        y = epochStorage(k, :) - mean(epochStorage(k, 1:250));
        psthMatrix(k,:) = binData(y, binRate, sampleRate);
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

% [datesD,DI2] = sort(datetime(epochStartTimeD)); % use this for nori package sorting 
spikeMatrix = spikeMatrix(DI2,:);


spikeMatrixUnbinned = spikeMatrixUnbinned(DI2,:);
% epochNumResponses = epochNumResponses(DI2);
psthMatrix = psthMatrix(DI2,:);
frameDwell = frameDwell(DI2);
barOrientation = barOrientation(DI2);
 barWidth = barWidth(DI2);
 apRadius = apRadius(DI2);
 bgClass = bgClass(DI2);
 seedList = seedList(DI2);
 noiseClass = noiseClass(DI2);
 barFrameDwell = barFrameDwell(DI2);
%  binaryIndex = find(contains(noiseClass,'binary'));
%  gaussianIndex = find(contains(noiseClass,'gaussian'));
 binaryIndex = contains(noiseClass,'binary');
gaussianIndex = contains(noiseClass,'gaussian');

 bgClassOrig = bgClass;
% %  
%  gaussianIndex=gaussianIndex'; 
%  binaryIndex = binaryIndex';
%  noiseClass = noiseClass';
%  bgClass = bgClass';
%  apRadius = apRadius';
%  barWidth = barWidth';
%  barOrientation = barOrientation';
%  barFrameDwell = barFrameDwell';
%  frameDwell = frameDwell';
%  seedList = seedList';
%% psth split\



% motionPSTH = psthMatrix(bgClass(gaussianIndex)==1,:);
% randomPSTH = psthMatrix(bgClass(gaussianIndex)==2,:);
% staticPSTH = psthMatrix(bgClass(gaussianIndex)==3,:);

motionPSTH = psthMatrix(bgClass(binaryIndex)==1,:);
randomPSTH = psthMatrix(bgClass(binaryIndex)==2,:);
staticPSTH = psthMatrix(bgClass(binaryIndex)==3,:);

% motionPSTH = psthMatrix(bgClass(barOrientation==330)==1,:);
% randomPSTH = psthMatrix(bgClass(barOrientation==330)==2,:);
% staticPSTH = psthMatrix(bgClass(barOrientation==330)==3,:);


% motionPSTH = psthMatrix(bgClass(binaryIndex==false)==1,:);
% randomPSTH = psthMatrix(bgClass(binaryIndex==false)==2,:);
% staticPSTH = psthMatrix(bgClass(binaryIndex==false)==3,:);


%repeated seed
% motionPSTH = psthMatrix(bgClass(seedList<=1)==1,:);
% randomPSTH = psthMatrix(bgClass(seedList<=1)==2,:);
% staticPSTH = psthMatrix(bgClass(seedList<=1)==3,:);



motionPSTH = mean(motionPSTH)';
randomPSTH = mean(randomPSTH)';
staticPSTH = mean(staticPSTH)';

msecX = linspace(0,10,10500);
plot(msecX,motionPSTH,'r','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 10 0 80])
% makeAxisStruct(gca,'seqPSTHSmooth')
figure
plot(msecX,randomPSTH,'b','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 10 0 80])
% makeAxisStruct(gca,'randPSTHSmooth')
figure
plot(msecX,staticPSTH,'k','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 10 0 80])
% makeAxisStruct(gca,'staticPSTHSmooth')

figure
plot(msecX,motionPSTH,'r','LineWidth',2)
hold on
plot(msecX,randomPSTH,'b','LineWidth',2)
plot(msecX,staticPSTH,'k','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 10 0 80])




figure
plot(msecX,motionPSTH,'r','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([5 10 0 80])
line([5 10], [mean(motionPSTH(5000:end,1)) mean(motionPSTH(5000:end,1))])
makeAxisStruct(gca,'seqPSTHBTzoomed')
figure
plot(msecX,randomPSTH,'b','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([5 10 0 80])
line([5 10], [mean(randomPSTH(5000:end,1)) mean(randomPSTH(5000:end,1))])
makeAxisStruct(gca,'randPSTHBTzoomed')
figure
plot(msecX,staticPSTH,'k','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([5 10 0 80])
line([5 10], [mean(staticPSTH(5000:end,1)) mean(staticPSTH(5000:end,1))])
makeAxisStruct(gca,'staticPSTHBTzoomed')

figure
plot(msecX,motionPSTH,'r','LineWidth',2)
hold on
plot(msecX,randomPSTH,'b','LineWidth',2)
plot(msecX,staticPSTH,'k','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 5 0 80])

title('all 3')
% save('psth3.mat','motionPSTH','randomPSTH','staticPSTH')

figure
plot(msecX,motionPSTH,'r','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 5 0 80])
line([0 5], [mean(motionPSTH(1:5000,1)) mean(motionPSTH(1:5000,1))])
makeAxisStruct(gca,'seqPSTHBTzoomedEARLY')
figure
plot(msecX,randomPSTH,'b','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 5 0 80])
line([0 5], [mean(randomPSTH(1:5000,1)) mean(randomPSTH(1:5000,1))])
makeAxisStruct(gca,'randPSTHBTzoomedEARLY')
figure
plot(msecX,staticPSTH,'k','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 5 0 80])
line([0 5], [mean(staticPSTH(1:5000,1)) mean(staticPSTH(1:5000,1))])
makeAxisStruct(gca,'staticPSTHBTzoomedEARLY')

figure
plot(msecX,motionPSTH,'r','LineWidth',2)
hold on
plot(msecX,randomPSTH,'b','LineWidth',2)
plot(msecX,staticPSTH,'k','LineWidth',2)
xlabel('Time (s)','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)
axis([0 5 0 80])
title('all 3')




%%

motionEpoch = epochStorage(bgClass(seedList<=1)==1,:);
randomEpoch = epochStorage(bgClass(seedList<=1)==2,:);
staticEpoch = epochStorage(bgClass(seedList<=1)==3,:);

figure
msecX2 = linspace(1,10500,105000);
plot(msecX2,motionEpoch(1,:),'r','LineWidth',2)
xlabel('Time (ms)','FontSize',24)
% ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'visible','off')
set(gca,'FontSize',20)
makeAxisStruct(gca,'seqSpikesSmooth')
figure
msecX2 = linspace(1,10500,105000);
plot(msecX2,randomEpoch(1,:),'b','LineWidth',2)
xlabel('Time (ms)','FontSize',24)
% ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'visible','off')
set(gca,'FontSize',20)
makeAxisStruct(gca,'randSpikesSmooth')
figure
msecX2 = linspace(1,10500,105000);
plot(msecX2,staticEpoch(1,:),'k','LineWidth',2)
xlabel('Time (ms)','FontSize',24)
% ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'visible','off')
set(gca,'FontSize',20)
makeAxisStruct(gca,'staticSpikesSmooth')
%% excel workaround

xlswrite('psthData.xls',motionPSTH,'A1')
xlswrite('psthData.xls',randomPSTH,'A2')
xlswrite('psthData.xls',staticPSTH,'A3')
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

gaussianOrBinary = 'gaussian';
% gaussianOrBinary = 'binary';
splitter2 = 90;
splitter3 = 350;
splitter4 = 1;
noiseSplitter = strcmp(gaussianOrBinary,noiseClass);

splitIndex2 = barOrientation == splitter2;
splitIndex3 = apRadius == splitter3;
% splitIndex3 = frameDwell == splitter3;
% splitIndex4 = barFrameDwell == splitter4;
splitIndex4 = seedList ~= splitter4;

dumbSplit = splitIndex2./splitIndex3./splitIndex4./noiseSplitter;
dumbSplit = dumbSplit==1;

% splitIndex3 = barWidth == splitter3;
% splitIndex4 = barFrameDwell == splitter4;bar
% splitIndex = splitIndex2 == splitIndex3(splitIndex3 ==0 ; 
% splitIndex = splitIndex == splitIndex4;
% bgClass = bgClass(splitIndex2>0);
% binaryIndex = dumbSplit==1;  
gaussianIndex = dumbSplit ==1;
% frameDwell = frameDwell(splitIndex2);

%% fake logicals
gaussianIndex(10:end) = false;

%% mean
plot(mean(lfRespStatic),'Color','k')
hold on
plot(mean(lfRespRandom),'Color','b')
plot(mean(lfRespSeq),'Color','r')
legend('static','random','sequential')
%% test only repeated sequence

if noiseFlag == 1 
noiseVars = struct();
noiseVars.type = 'gaussian';
noiseVars.contrast = 0.3333;

timings = [250,10000,250]; % AUTOMATE THIS LATER
frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage(gaussianIndex,:), 10000, binRate);
% frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage, 10000, 10000);
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,2,seedList(gaussianIndex),frameDwell(gaussianIndex));
response1 = psthMatrix(gaussianIndex,:);
else
    noiseVars = struct();
noiseVars.type = 'binary';
noiseVars.contrast = 0.3333;

timings = [250,10000,250]; % AUTOMATE THIS LATER
frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage(binaryIndex,:), 10000, binRate);
% frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage, 10000, 10000);
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,2,seedList(binaryIndex),frameDwell(binaryIndex));
response1 = psthMatrix(binaryIndex,:);
end

%      
     
 frameValues = frameValuesAll;
%  findCondition = size(response1,1);
%  numRuns = findCondition/3;
%  seqRuns = 1:3:numRuns*3;
%  randRuns = 2:3:numRuns*3;
%  staticRuns= 3:3:numRuns*3;


% ftest = ifft(fft(response1(3,:)) .* conj(fft(frameValues(1,:))))
% plot(ftest)
% 
% rSeedGenerator = ifft(fft(linearFilterAvg) .* fft(frameValues(3,:)));
% rSeedY = outputNonlinearity(nlParams,rSeedGenerator);
% plot(response1(3,:));
% hold on 
% plot(rSeedY);
    

 %% Get stim frames -- #5 (standard, using Mike's frame time functions)
count = 0;
noiseFlag = 1; %1 for gaussian
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
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,2,seedList(gaussianIndex),frameDwell(gaussianIndex));
else
    noiseVars = struct();
noiseVars.type = 'binary';
noiseVars.contrast = 0.3333;

timings = [250,10000,250]; % AUTOMATE THIS LATER
frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage(binaryIndex,:), 10000, binRate);
% frames = manookinlab.ovation.getFrameTimesFromMonitor(monitorStorage, 10000, 10000);
frameValuesAll = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,2,seedList(binaryIndex),frameDwell(binaryIndex));
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

lfRespSeq(countSeq,:) = lfResp;
stimOutSeq(countSeq,:) = lfStim;

    elseif bgClass(h) == 2
        countR = countR+1;
lfResp = response1(h,:);
lfStim = frameValues(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterRandom = (linearFilterRandom*(countR-1) + lf)/countR;
sortNLRand(countR) = h;

lfRespRandom(countR,:) = lfResp;
stimOutRandom(countR,:) = lfStim;

    elseif bgClass(h) == 3
        countStatic = countStatic + 1;
lfResp = response1(h,:);
lfStim = frameValues(h,:);
lf = real(ifft( fft(lfResp(:)') .* conj(fft(lfStim(:)')) ));

linearFilterStatic = (linearFilterStatic*(countStatic-1) + lf)/countStatic;
sortNLStatic(countStatic) = h;

lfRespStatic(countStatic,:) = lfResp;
stimOutStatic(countStatic,:) = lfStim;

    end
        
end

% save for decomp

save('decompstuff.mat','lfRespSeq','lfRespRandom','lfRespStatic','stimOutSeq','stimOutRandom','stimOutStatic');

%normalize to scale the filters 
linearFilterSeq = linearFilterSeq/norm(linearFilterSeq);
linearFilterRandom = linearFilterRandom/norm(linearFilterRandom);
linearFilterStatic = linearFilterStatic/norm(linearFilterStatic);

linearFilterAvg = (linearFilterSeq + linearFilterRandom + linearFilterStatic)/3;

%fits for filters using Obsidian functions 
     seqParams = models.ln.fitLinearFilterParams(linearFilterSeq,binRate);
     filterSeq = linearFilterFunction(seqParams,(1:plotLngth)/binRate);
     figure(7)
     plot((1:plotLngth)/binRate, filterSeq,'LineWidth',2,'Color','r')
     hold on
     
     randomParams = models.ln.fitLinearFilterParams(linearFilterRandom,binRate);
     filterRand = linearFilterFunction(randomParams,(1:plotLngth)/binRate);
     plot((1:plotLngth)/binRate, filterRand,'LineWidth',2,'Color','b')
     
     staticParams = models.ln.fitLinearFilterParams(linearFilterStatic,binRate);
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
line((1:plotLngth)/binRate, linearFilterAvg(1:plotLngth),'LineStyle',':','LineWidth',1,'Color','magenta')
set(gca,'xdir','reverse')
title('linear filters for each surround condition')
legend('Sequential','Random','Static','Average','Location','northwest')
xlabel('time(s)')
ylabel('weight')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
% export_fig 'noiseNL.pdf'

filterX = (1:plotLngth)/binRate;
filterX=filterX';
lFilter=[];
lFilter(1,:)=linearFilterSeq(1:plotLngth);
lFilter(2,:)=linearFilterRandom(1:plotLngth);
lFilter(3,:)=linearFilterStatic(1:plotLngth);

lFilter=lFilter';
save('PoffOuttFilter','lFilter','filterX')

     
 
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


% pred = ifft(fft(seqFV(p,:)) .* fft(linearFilterSeq(:)'));
pred = ifft(fft(seqFV(p,:)) .* fft(linearFilterAvg(:)'));
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
 
opts = statset('nlinfit');
opts.MaxFunEvals = 1e5;
opts.MaxIter = 1e5;
 
g = nlinfit(huhx,huhy,@outputNonlinearity,[max(yBin)*3 0.2 -1.5 min(yBin)],opts);
% nonLinXAx1= linspace(min(huhx),max(huhx),1000);
nonLinXAx1 = linspace(-1,1,1000);
nonLinXAx2 = linspace(-1,1,100);
% nonLinXAx3= linspace(-1,1,1000);
nonLinXAx3= linspace(min(huhx),max(huhx),1000);

 
  

 figure(8);
 clf
 plot(huhx/max(huhx),huhy,'.','Color','r')

 
%  plot(nonLinXAx3,huhy*10000,'.','Color','r')
 hold on
 plot(nonLinXAx3/max(nonLinXAx3),outputNonlinearity(g,nonLinXAx3),'Color','r','LineWidth',1.5)
axis([-1 1 0 max(huhy)+2])

  nonLinFit = outputNonlinearity(g,nonLinXAx3);
  


%  predStore(:,stimExtent(501:end))
[xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(:)),respStore(:,stimExtent(:)),nonlinearityBins);
 nlParams = fitNonlinearityParams(xfBin, yfBin);

 nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),1000);
 figure(9)
 figure(98)
 plot(nlX,outputNonlinearity(nlParams,nlX),'Color','r','LineWidth',2)
 hold on
 nlParamsSeq = nlParams;
 outY=[];
 outX=[];
 outY(1,:)=outputNonlinearity(nlParams,nlX);
 outX(1,:)=nlX;
 
 
 
 
 binnedX(:,1) = xfBin;
 binnedY(:,1) = yfBin; 
 
 XforModel(:,1) = nlX;
 
 
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
% pred = ifft(fft(randFV(p,:)) .* fft(linearFilterRandom(:)'));
pred = ifft(fft(randFV(p,:)) .* fft(linearFilterAvg(:)'));
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
 


g = nlinfit(huhx,huhy,@outputNonlinearity,[max(yBin)*3 0.2 -1.5 min(yBin)],opts);
% nonLinXAx1= linspace(min(huhx),max(huhx),1000);
nonLinXAx3= linspace(min(huhx),max(huhx),1000);






 figure(8)
 plot(huhx/max(huhx),huhy,'.','Color','b') %huhx*10
 

 

 
%  plot(nonLinXAx3,huhy*10000,'.','Color','b')
 hold on
 plot(nonLinXAx3/max(nonLinXAx3),outputNonlinearity(g,nonLinXAx3),'Color','b','LineWidth',1.5)
 
   nonLinFit = outputNonlinearity(g,nonLinXAx3);
  

 
%  axis([-1 1 0 max(huhy*10000)+2])
 [xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(:)),respStore(:,stimExtent(:)),nonlinearityBins);
%  [xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(1:9770)),respStore(:,stimExtent(1:9770)),nonlinearityBins);
 nlParams = fitNonlinearityParams(xfBin, yfBin);
 nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),1000);
 figure(9)
 figure(98)
 plot(nlX,outputNonlinearity(nlParams,nlX),'Color','b','LineWidth',2)
 outY(2,:)=outputNonlinearity(nlParams,nlX);
 outX(2,:)=nlX;
 hold on
 nlParamsRandom = nlParams;
 
 binnedX(:,2) = xfBin;
 binnedY(:,2) = yfBin;
 
 XforModel(:,2) = nlX;

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
% pred = ifft(fft(staticFV(p,:)) .* fft(linearFilterStatic(:)'));
pred = ifft(fft(staticFV(p,:)) .* fft(linearFilterAvg(:)'));
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
 
g = nlinfit(huhx,huhy,@outputNonlinearity,[max(yBin)*3 0.2 -1.5 min(yBin)],opts);
nonLinXAx3= linspace(min(huhx),max(huhx),1000);


     
 

 
 figure(8)
 plot(huhx/max(huhx),huhy,'.','Color','k')
 


 hold on
 plot(nonLinXAx3/max(nonLinXAx3),outputNonlinearity(g,nonLinXAx3),'Color','k','LineWidth',1.5)
  axis([-1 1 0 max(huhy)+2])
  
  nonLinFit = outputNonlinearity(g,nonLinXAx3);
  

 
%   axis([0 400 0 4*10^-3])
 hold on
 legend('sequential NL','sequential fit','random NL','random fit','static NL','static fit','location','northwest')
 title('NL for each surround condition')
 xlabel('input','FontSize',18)
 ylabel('spikes/s','FontSize',18)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
% export_fig 'noiseNL.pdf' -append

[xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(:)),respStore(:,stimExtent(:)),nonlinearityBins);
%  [xfBin,yfBin] = binNonlinearity(pred2Store(:,stimExtent(1:9770)),respStore(:,stimExtent(1:9770)),nonlinearityBins);
 nlParams = fitNonlinearityParams(xfBin, yfBin);
 nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),1000);
 figure(9)
 figure(98)
 plot(nlX,outputNonlinearity(nlParams,nlX),'Color','k','LineWidth',2)
  legend('sequential','random','static','location','northwest')
  xlabel('input')
  ylabel('spikes/s')
  
   outY(3,:)=outputNonlinearity(nlParams,nlX);
 outX(3,:)=nlX;
%  
 outY=outY';
 outX=outX';
 
 save('OffPout','outY','outX')
  
   binnedX(:,3) = xfBin;
   binnedY(:,3) = yfBin; 
 
   starterParams = [.2,-1.5];
% nlX2 = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),100);
 
mParams = fitMultiVarParams(binnedX(:,[1 2]),binnedY(:,[1 2]),1,starterParams);
% outNL = multiVarNL(mParams,XforModel);
outNL = multiHZNL(mParams,binnedX);

% nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),size(outNL,1))

figure
plot(binnedX(:,1),outNL(:,1),'Color',[.5 0 0],'LineStyle','-.','LineWidth',1.5)
hold on
plot(binnedX(:,2),outNL(:,2),'Color',[0 0 .5],'LineStyle','-.','LineWidth',1.5)
plot(binnedX(:,1),binnedY(:,1),'r','LineWidth',.5)
plot(binnedX(:,2),binnedY(:,2),'k','LineWidth',.5)
% plot(binnedX(:,2),binnedY(:,2),'b','LineWidth',.5)
xlabel('input')
ylabel('output')
title('x shift')

% legend('fit seq','fit rand','data seq','data rand')
legend('fit seq','fit static','data seq','data static')

% R2(binnedY(:,1),outNL(:,1))
% 
% R2(binnedY(:,2),outNL(:,2))
% 
% [rsquaredOut,rootMSE]=rsquare(binnedY(:,1),outNL(:,1))
% 
% [rsquaredOut,rootMSE]=rsquare(binnedY(:,2),outNL(:,2))

MSE1 = immse(binnedY(:,1),outNL(:,1))
MSE2 = immse(binnedY(:,2),outNL(:,2))

%save! for fit figure
youtSeq = binnedY(:,1);
youtRand = binnedY(:,2);
youtStatic = binnedY(:,3);
xoutSeq = binnedX(:,1);
xoutRand = binnedX(:,2);
xoutStatic = binnedX(:,3);
xoutSeqFitHZ = outNL(:,1);
xoutRandFitHZ = outNL(:,2);
 
mParams = fitMultiVarParams(binnedX(:,[1 2]),binnedY(:,[1 2]),0,starterParams);
% outNL = multiGainNL(mParams,XforModel);
outNL = multiGainNL(mParams,binnedX);


xoutSeqFitGain = outNL(:,1);
xoutRandFitGain = outNL(:,2);

save('NLandFit.mat','youtSeq','youtRand','youtStatic','xoutSeq','xoutRand','xoutStatic','xoutSeqFitHZ','xoutRandFitHZ','xoutSeqFitGain','xoutRandFitGain');


% nlX = linspace(-max(abs(xfBin(:))),max(abs(xfBin(:))),size(outNL,1))

figure
plot(binnedX(:,1),outNL(:,1),'Color',[.5 0 0],'LineStyle','-.','LineWidth',1.5)
hold on
plot(binnedX(:,2),outNL(:,2),'Color',[0 0 .5],'LineStyle','-.','LineWidth',1.5)
plot(binnedX(:,1),binnedY(:,1),'r','LineWidth',.5)
% plot(binnedX(:,2),binnedY(:,2),'b','LineWidth',.5)
plot(binnedX(:,2),binnedY(:,2),'k','LineWidth',.5)
xlabel('input')
ylabel('output')
title('gain change')

% legend('fit seq','fit rand','data seq','data rand')
legend('fit seq','fit static','data seq','data static')

% R2(binnedY(:,1),outNL(:,1))

% R2(binnedY(:,2),outNL(:,2))

% [rsquaredOut,rootMSE]=rsquare(binnedY(:,1),outNL(:,1))

% [rsquaredOut,rootMSE]=rsquare(binnedY(:,2),outNL(:,2))

MSE1 = immse(binnedY(:,1),outNL(:,1))
MSE2 = immse(binnedY(:,2),outNL(:,2))

  

% if saveFlag == 1 && strcmp(cellType,'OFF Parasol')
%     MNOFFParasol.Exp(expNum).Cell(OFFParasolC).NL.Static = outputNonlinearity(nlParams,nlX);
% elseif saveFlag == 1 && strcmp(cellType,'ON Parasol')
%     MNONParasol.Exp(expNum).Cell(ONParasolC).NL.Static = outputNonlinearity(nlParams,nlX);
% elseif saveFlag ==1 && strcmp(cellType,'ON Smooth')
%     MNONSmooth.Exp(expNum).Cell(ONSmoothC).NL.Static = outputNonlinearity(nlParams,nlX);
% elseif saveFlag ==1 && strcmp(cellType,'OFF Smooth')
%     MNOFFSmooth.Exp(expNum).Cell(OFFSmoothC).NL.Static = outputNonlinearity(nlParams,nlX);
% end



%  x=-300:400;
% nlfun = @(p,x)(p(1)*normcdf(p(2)*x+p(3),0,1));
% p = nlinfit(huhx,huhy,nlfun,[5 50 0]);
% 
% 
% 
% plot(xBin,yBin,'.');
% plot(x,nlfun(p,x));
% hold off;
figure(98)
plot(binnedX(:,1),binnedY(:,1),'.','Color','r')
hold on
plot(binnedX(:,2),binnedY(:,2),'.','Color','b')
plot(binnedX(:,3),binnedY(:,3),'.','Color','k')
% title('binned NLs')
legend('Motion','Random','Static','location','NorthWest')
xlabel('Input','FontSize',24)
ylabel('Spike Rate (Hz)','FontSize',24)
set(gca,'box','off') 
set(gca,'FontSize',20)



figure(99)
plot(binnedX(:,1),binnedY(:,1),'Color','r')
hold on
plot(binnedX(:,2),binnedY(:,2),'Color','b')
plot(binnedX(:,3),binnedY(:,3),'Color','k')
title('binned NLs')
legend('Motion','Random','Static','location','NorthWest')

%% final save
cd('/Users/toddappleby/Documents/Data/Clarinet Exports/SavedData')
save('MotionandNoise.mat','MNOFFParasol')
save('MotionandNoise.mat','-append','MNONParasol')
save('MotionandNoise.mat','-append','MNONSmooth')
save('MotionandNoise.mat','-append','MNOFFSmooth')

%% Motion CS
dataType =1;
count = 0;
desiredSTD = 6;
clear epochStorage 
clear bgClass
clear centerClass
clear contrastState
count2 = 0;
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
%    recordingTechnique = epochs(i).meta.recordingTechnique;
   egLabel = epochs(i).meta.epochGroupLabel;
   
   if strcmp(displayName,'Motion Center Surround')
       count2 = count2+1;
       
       egLabel = epochs(i).meta.epochGroupLabel;
       surroundO(count2) = epochs(i).meta.surroundBarOrientation;
       centerO(count2) = epochs(i).meta.centerBarOrientation;
        width = epochs(i).meta.surroundBarWidth;
        cWidth = epochs(i).meta.centerBarWidth;
          surroundOrientation = epochs(i).meta.surroundBarOrientation;
          centerBarOrientation = epochs(i).meta.centerBarOrientation;
          surroundBarDwell = epochs(i).meta.surroundBarFrameDwell;
          centerBarDwell = epochs(i).meta.centerFrameDwell;
          motionTime = epochs(i).meta.motionTime;
          numBars = epochs(i).meta.numberOfCenterBars;
        
        if  surroundOrientation == 90 && centerBarOrientation == 90 && strcmp(egLabel,'Control')
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
[meanCollection, meanKey] = motionCS(epochs,bgClass,centerClass,contrastState,epochStorage,protocolExample,saveGraph,desiredSTD,dataType);



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
%         psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
psthMatrix(k,:) = psth(spikeMatrix(k,:),10,sampleRate,1);
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

function [meanCollection, meanKey] = motionCS(epochs,bgClass,centerClass,contrastState,epochStorage,protocolExample,saveGraph,STD1,dataType)
% dataType = 1;
sampleRate = 10000;
binRate = 1000;
preTime = epochs(protocolExample).meta.preTime;
stimTime = epochs(protocolExample).meta.stimTime;
stimOrig = stimTime;
motionTime = epochs(protocolExample).meta.motionTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;
desiredSTD = STD1;
% spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
% psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
if dataType == 1
for k = 1:size(epochStorage,1)
    bgClass(k,1);
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else
              
%         spikes = binSpikeCount(spikes, binRate, sampleRate);
%         spikes = spikes/max(spikes);
        spikes = spikes';
        spikeMatrix(k,:) = spikes; 
        
        psthMatrix(k,:) = psth(spikeMatrix(k,:),20,sampleRate,1);
%         spikeMatrix(k,:) = spikes;
%         psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end
else
for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
end
end

size(spikeMatrix)
size(psthMatrix)
figure(90)
plot(spikeMatrix(12,:))

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
preTime = preTime * 10;
motionTime=motionTime*10;

motionStart = 1000; % use these to extend/shorten spike count window
motionEnd = 2500;
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

a=sum(psthMatrix(posContrast(posSeq),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
b=sum(psthMatrix(posContrast(posRand),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
% a=a/.1;
% b=b/.1;


PCS = mean(a);
PCSerror = sem(a);

PCR = mean(b);
PCRerror = sem(b);

c=sum(psthMatrix(negContrast(negSeq),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
d=sum(psthMatrix(negContrast(negRand),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
% c=c/.1;
% d=d/.1;



NCS = mean(c);
NCSerror=sem(c);

NCR = mean(d);
NCRerror = sem(d);

e=sum(psthMatrix(posContrast(centerSeq),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
f=sum(psthMatrix(posContrast(centerRand),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
g=sum(psthMatrix(posContrast(center180),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);

% e=e/.1;
% f=f/.1;
% g=g/.1;


CPCS= mean(e);
CPCSerror = sem(e);
CPCR = mean(f);
CPCRerror = sem(f);
CPC180 = mean(g);
CPC180error = sem(g);

h=sum(psthMatrix(posContrast(centerRSeq),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
i=sum(psthMatrix(posContrast(centerRRand),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
j=sum(psthMatrix(posContrast(centerR180),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);

h=h;
i=i;
j=j;

rCPCS= mean(h);
rCPCSerror = sem(h);
rCPCR = mean(i);
rCPCRerror = sem(i);
rCPC180 = mean(j);
rCPC180error = sem(j);

k=sum(psthMatrix(negContrast(centerSeqNeg),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
l=sum(psthMatrix(negContrast(centerRandNeg),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
m=sum(psthMatrix(negContrast(center180Neg),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);

k=k;
l=l;
m=m;

CNCS= mean(k);
CNCSerror = sem(k);
CNCR = mean(l);
CNCRerror = sem(l);
CNC180 = mean(m);
CNC180error = sem(m);

n=sum(psthMatrix(negContrast(centerRSeqNeg),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
o=sum(psthMatrix(negContrast(centerRRandNeg),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);
p=sum(psthMatrix(negContrast(centerR180Neg),((preTime)+(motionTime)-motionStart):((preTime)+(motionTime)+(motionEnd))),2);

n=n;
o=o;
p=p;

rCNCS= mean(n);
rCNCSerror = sem(n);
rCNCR = mean(o);
rCNCRerror = sem(o);
rCNC180 = mean(p);
rCNC180error = sem(p);

meanCollection = [PCS PCR NCS NCR CPCS CPCR CPC180 rCPCS rCPCR rCPC180 CNCS CNCR CNC180 rCNCS rCNCR rCNC180];
meanKey = ['PCS ', 'PCR ', 'NCS ', 'NCR ', 'CPCS ', 'CPCR ', 'CPC180 ', 'rCPCS ', 'rCPCR ', 'rCPC180 ', 'CNCS ', 'CNCR ', 'CNC180 ', 'rCNCS ', 'rCNCR ', 'rCNC180 '];

barGraphArray = [CPCS CPCR CPC180 rCPCS rCPCR rCPC180 CNCS CNCR CNC180 rCNCS rCNCR rCNC180]';

barGraphError = [CPCSerror CPCRerror CPC180error rCPCSerror rCPCRerror rCPC180error CNCSerror CNCRerror CNC180error rCNCSerror rCNCRerror rCNC180error]';

save('barGraphStuff.mat','barGraphArray','barGraphError')

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
% posContrastSeq(:,1) = mean(psthMatrix(posContrast(centerSeq),:));
% posContrastSeq(:,2) = mean(psthMatrix(posContrast(centerRand),:));
% posContrastSeq(:,3) = mean(psthMatrix(posContrast(center180),:));
% save(['E:\Data Analysis_2020\weeklymeeting_0522\posContrastSeq.mat'],'posContrastSeq');


subplot(2,3,4)
plot(mean(psthMatrix(posContrast(centerRSeq),:)),'Color','r','LineWidth',1)
hold on
plot(mean(psthMatrix(posContrast(centerRRand),:)),'Color','b','LineWidth',1)
plot(mean(psthMatrix(posContrast(centerR180),:)),'Color','k','LineWidth',1)
legend('seq','rand','180','location','northwest')
title('all positive center, random surround')
% posContrastRand(:,1) = mean(psthMatrix(posContrast(centerRSeq),:));
% posContrastRand(:,2) = mean(psthMatrix(posContrast(centerRRand),:));
% posContrastRand(:,3) = mean(psthMatrix(posContrast(centerR180),:));
% save(['E:\Data Analysis_2020\weeklymeeting_0522\posContrastRand.mat'],'posContrastRand');


subplot(2,3,5)
plot(mean(psthMatrix(negContrast(centerSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(center180Neg),:)),'Color','k')
legend('seq','rand','180','location','northwest')
title('all negative center, sequential surround')
% negContrastSeq(:,1) = mean(psthMatrix(negContrast(centerSeqNeg),:));
% negContrastSeq(:,2) = mean(psthMatrix(negContrast(centerRandNeg),:));
% negContrastSeq(:,3) = mean(psthMatrix(negContrast(center180Neg),:));
% save(['E:\Data Analysis_2020\weeklymeeting_0522\negContrastSeq.mat'],'negContrastSeq');

subplot(2,3,6)
plot(mean(psthMatrix(negContrast(centerRSeqNeg),:)),'Color','r')
hold on
plot(mean(psthMatrix(negContrast(centerRRandNeg),:)),'Color','b')
plot(mean(psthMatrix(negContrast(centerR180Neg),:)),'Color','k')
legend('seq','rand','180','location','northwest')
title('all negative center, random surround')
% negContrastRand(:,1) = mean(psthMatrix(negContrast(centerRSeqNeg),:));
% negContrastRand(:,2) = mean(psthMatrix(negContrast(centerRRandNeg),:));
% negContrastRand(:,3) = mean(psthMatrix(negContrast(centerR180Neg),:));
% save(['E:\Data Analysis_2020\weeklymeeting_0522\negContrastRand.mat'],'negContrastRand');

figure(70)

subplot(1,2,1)

plot(mean(psthMatrix(negContrast(centerSeqNeg),:))','Color','r')

negSeqSeq = mean(psthMatrix(negContrast(centerSeqNeg),:))';

hold on

plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')

negRandSeq = mean(psthMatrix(negContrast(centerRandNeg),:))';

plot(mean(psthMatrix(negContrast(center180Neg),:)),'Color','k')

neg180Seq = mean(psthMatrix(negContrast(center180Neg),:))';

ylabel('Spike/sec','FontSize',14)

xlabel('Time (ms)','Fontsize',14)

legend('Motion','Random','Motion-Reversed','location','northwest')

title('3 Center Conditions w/ Surround Motion','FontSize',12)

 

 

subplot(1,2,2)

plot(mean(psthMatrix(negContrast(centerRSeqNeg),:)),'Color','r')

negSeqRand = mean(psthMatrix(negContrast(centerRSeqNeg),:))';

hold on

plot(mean(psthMatrix(negContrast(centerRRandNeg),:)),'Color','b')

negRandRand = mean(psthMatrix(negContrast(centerRRandNeg),:))';

plot(mean(psthMatrix(negContrast(centerR180Neg),:)),'Color','k')

neg180Rand = mean(psthMatrix(negContrast(centerR180Neg),:))';




ylabel('Spike/sec','FontSize',14)

xlabel('Time (ms)','Fontsize',14)

legend('Motion','Random','Motion-Reversed','location','northwest')

title('3 Center Conditions w/ Non-motion Surround','FontSize',12)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')


% figure(50)
% subplot(1,2,1)
% plot(mean(psthMatrix(negContrast(centerSeqNeg),:)),'Color','r')
% hold on
% plot(mean(psthMatrix(negContrast(centerRandNeg),:)),'Color','b')
% plot(mean(psthMatrix(negContrast(center180Neg),:)),'Color','k')
% legend('seq','rand','180','location','northwest')
% title('all negative center, sequential surround')
% 
% subplot(1,2,2)
% plot(mean(psthMatrix(negContrast(centerRSeqNeg),:)),'Color','r')
% hold on
% plot(mean(psthMatrix(negContrast(centerRRandNeg),:)),'Color','b')
% plot(mean(psthMatrix(negContrast(centerR180Neg),:)),'Color','k')
% legend('seq','rand','180','location','northwest')
% title('all negative center, random surround')

figure(51)
subplot(1,2,1)

plot(mean(psthMatrix(posContrast(centerSeq),:)),'Color','r','LineWidth',1)
% -mean(mean(psthMatrix(posContrast(centerSeq),1:5000)))
posSeqSeq = mean(psthMatrix(posContrast(centerSeq),:))';

hold on

plot(mean(psthMatrix(posContrast(centerRand),:)),'Color','b','LineWidth',1)
% -mean(mean(psthMatrix(posContrast(centerRand),1:5000)))
posRandSeq = mean(psthMatrix(posContrast(centerRand),:))';

plot(mean(psthMatrix(posContrast(center180),:)),'Color','k','LineWidth',1)
% -mean(mean(psthMatrix(posContrast(center180),1:5000)))
pos180Seq= mean(psthMatrix(posContrast(center180),:))';

ylabel('Spike/sec','FontSize',14)

xlabel('Time (ms)','Fontsize',14)

legend('Motion','Random','Motion-Reversed','location','northwest')

title('3 Center Conditions w/ Surround Motion','FontSize',12)

subplot(1,2,2)
plot(mean(psthMatrix(posContrast(centerRSeq),:)),'Color','r','LineWidth',1)
% -mean(mean(psthMatrix(posContrast(centerRSeq),1:5000)))
posSeqRand = mean(psthMatrix(posContrast(centerRSeq),:)');

hold on
plot(mean(psthMatrix(posContrast(centerRRand),:)),'Color','b','LineWidth',1)
% -mean(mean(psthMatrix(posContrast(centerRRand),1:5000)))
posRandRand = mean(psthMatrix(posContrast(centerRRand),:))';

plot(mean(psthMatrix(posContrast(centerR180),:)),'Color','k','LineWidth',1)
% -mean(mean(psthMatrix(posContrast(centerR180),1:5000)))
pos180Rand = mean(psthMatrix(posContrast(centerR180),:))';

ylabel('Spike/sec','FontSize',14)

xlabel('Time (ms)','Fontsize',14)

legend('Motion','Random','Motion-Reversed','location','northwest')

title('3 Center Conditions w/ Non-motion Surround','FontSize',12)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')

save('mcsPSTH.mat','posSeqSeq','posRandSeq','pos180Seq','posSeqRand','posRandRand','pos180Rand','negSeqSeq','negRandSeq','neg180Seq','negSeqRand','negRandRand','neg180Rand')

scaleX = linspace(1,2000,20000);
scaleX = scaleX';
save('scalex.mat','scaleX')


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