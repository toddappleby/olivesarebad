%% CHOOSE CELL TO LOAD
exportDirectory = 'C:\Users\reals\Documents\PhD 2021\ClarinetExports\'; %new PC 
experimentDate = '2023_0607';
cellNum = 'Ac3';
cd(strcat(exportDirectory,experimentDate)) 
%20210910Fc1
%% LOAD AND GRAB CRITICAL ELEMENTS

cd(strcat(exportDirectory,experimentDate))
load(strcat(experimentDate(1:4),experimentDate(6:9),cellNum,'_FT.mat'))
frameTs = epochs;
load(strcat(experimentDate(1:4),experimentDate(6:9),cellNum,'.mat'))
cellData = epochs;

%get list of protocols run, display one of each.
for z = 1:size(epochs,2)
   list(z) = string(epochs(z).meta.displayName);
end
uniqueProtocols = unique(list'); %no idea why i didn't do this to begin with, at least for how I ended up using the list

% define static parameters
sampleRate = 10000;
frameRate = 60;
binRate = 1000;

% grab frame timings. Faster if kept in cell array (cell2mat is slow)
frameTimings= [];
dNameLogical = strcmp(string(arrayfun(@(x)(getfield(x(1).meta,'displayName')), cellData, 'UniformOutput', false)),'Motion And Noise');
dNameLogicalFrames = strcmp(string(arrayfun(@(x)(getfield(x(1).meta,'displayName')), frameTs, 'UniformOutput', false)),'Motion And Noise'); %2019 thru October's files (at least) don't match between cellData and frameTs.  Haven't investigated why
frameTimings = {frameTs(dNameLogicalFrames).epoch};
frameTimings = cell2mat(frameTimings'); %needs to be double for use later


% get raw spiking data for each epoch
rawData= [];
rawData = {cellData(dNameLogical).epoch}; 

% get metadata for each epoch (a longer process)
metaData=[];
metaData = cellData(dNameLogical); %grab only protocol data

    % critical surround state info for motion and noise protocol
    bgClass= string(arrayfun(@(x)(getfield(x(1).meta,'backgroundClass')), metaData, 'UniformOutput', false));
    bgClass(strcmp(bgClass,"sequential"))="1"; bgClass(strcmp(bgClass,"random"))="2"; bgClass(strcmp(bgClass,"stationary"))="3"; % assign surround states a num for ease of use
    bgClass=double(bgClass);


    %individual parameters that vary and we care about.  could make more
    %flexible ... but WHY?
    seedList = double(string(arrayfun(@(x)(getfield(x(1).meta,'seed')), metaData, 'UniformOutput', false)));
    
    noiseClass = string(arrayfun(@(x)(getfield(x(1).meta,'noiseClass')), metaData, 'UniformOutput', false));
    apRadius = double(string(arrayfun(@(x)(getfield(x(1).meta,'apertureRadius')), metaData, 'UniformOutput', false)));
    barOrientation = double(string(arrayfun(@(x)(getfield(x(1).meta,'barOrientation')), metaData, 'UniformOutput', false)));
    epochLabel = string((string(arrayfun(@(x)(getfield(x(1).meta,'epochGroupLabel')), metaData, 'UniformOutput', false))));
    epochTime = arrayfun(@(x)(getfield(x(1).meta,'epochTime')), metaData, 'UniformOutput', false);
    
    if isfield(metaData(1).meta,'frameDwell') % old data frameDwell wasn't changeable parameter
    frameDwell = double(string(arrayfun(@(x)(getfield(x(1).meta,'frameDwell')), metaData, 'UniformOutput', false)));
    else
        frameDwell(1:length(metaData))=1;
    end
    
    %params that likely won't change:
    preTime = metaData(1).meta.preTime; %only need once because unlikely to change btwn epochs
    stimTime = metaData(1).meta.stimTime;
    tailTime = metaData(1).meta.tailTime;

        % convert time to datapoints
            prePts = preTime*10; % msec to pts
            stimPts = stimTime*10;
            tailpts = tailTime*10;
            
        %sort by time
        timeArray = [epochTime{1:end}];
        [~,timeSorter] = sort(timeArray);
        apRadius=apRadius(timeSorter);
        barOrientation=barOrientation(timeSorter);
        epochLabel=epochLabel(timeSorter);
        noiseClass = noiseClass(timeSorter);
        seedList = seedList(timeSorter);
        bgClass = bgClass(timeSorter);
        frameTimings = frameTimings(timeSorter,:);
        rawData = rawData(timeSorter); % cell array
        
        
            
%% Grab Spikes, 

desiredSTD = 4; %standard deviation used to detect spikes

for k = 1:length(rawData)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(rawData{k}, [prePts stimPts], ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else             
        spikeMatrixUnbinned(k,:) = spikes;
        spikes = binSpikeCount(spikes, binRate, sampleRate);
        spikes = spikes';
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:)*binRate,6+2/3,binRate,1);
          end
end

[uniques, ~, classLabel] = unique(noiseClass);
uniques(mode(classLabel))


'radius'
mode(apRadius)
'orientation'
mode(barOrientation)
'framedwell'
mode(frameDwell)
%% spike times for auto correlation ! 

allTimes = [];
onlyStatic = find(bgClass==3);


for t = onlyStatic
    
    spikeTimes = find(spikeMatrix(t,:)>0);
    isiTmp = diff(spikeTimes)';
    allTimes = [allTimes; isiTmp];
end

[y,x]= getSpikeAutocorrelation(allTimes); 

figure
plot(y)



%% Select Parameter Combinations As Desired

binaryIndex = contains(noiseClass,'binary')';
gaussianIndex = contains(noiseClass,'gaussian')';
% 
gaussianOrBinary = 'gaussian'; %choose your fighter
% gaussianOrBinary = 'binary';

egLabelCompare = "Control";
% egLabelCompare = "AltCenter";
% egLabelCompare = "motion and noise";


noiseSplitter = strcmp(gaussianOrBinary,noiseClass);
egSplitter = strcmp(egLabelCompare,epochLabel);
splitter2 = 90;
splitter3 = 350;
splitter4 = 1;
splitter5 = 1;

    %typical sorters
    splitIndex2 = barOrientation == splitter2;
    splitIndex3 = apRadius == splitter3;
    splitIndex4 = frameDwell == splitter4;
    splitIndex5 = seedList ~= splitter5; %change this to true for repeated seeds (   ) 

dumbSplit = splitIndex2./splitIndex3./splitIndex4./splitIndex5./noiseSplitter./egSplitter;
dumbSplit = dumbSplit==1;

if strcmp(gaussianOrBinary,'gaussian')
gaussianIndex = dumbSplit ==1;
usableIndex = gaussianIndex;
else
binaryIndex = dumbSplit==1; 
usableIndex = binaryIndex;
end
%% first go 
response = psthMatrix(usableIndex,:);
responseRaw = spikeMatrixUnbinned(usableIndex,:);
responseSpikes = spikeMatrix(usableIndex,:);
bgClassX = bgClass(usableIndex);

% seqMean = mean(response(bgClassX==1,:));
% randMean = mean(response(bgClassX==2,:));
% staticMean = mean(response(bgClassX==3,:));


seqMean = mean(response(bgClassX==1,:),1);
randMean = mean(response(bgClassX==2,:),1);
staticMean = mean(response(bgClassX==3,:),1);

seq2Mean = mean(responseRaw(bgClassX==1,:),1);
rand2Mean = mean(responseRaw(bgClassX==2,:),1);
static2Mean = mean(responseRaw(bgClassX==3,:),1);

figure(15)
plot(seqMean,'r')
hold on
plot(randMean,'b')
hold on
plot(staticMean,'k')
legend('sequential','random','static')

figure(17)
plot(seq2Mean,'r')
hold on
plot(rand2Mean,'b')
hold on
plot(static2Mean,'k')
legend('sequential','random','static')

%generate psth
figure(16)
subplot(3,1,1)
    plot(mean(response(bgClassX==3,:)),'k')
subplot(3,1,2)
    plot(mean(response(bgClassX==2,:)),'b')
subplot(3,1,3)
    plot(mean(response(bgClassX==1,:)),'r')

    figure
    plot(mean(response(bgClassX==3,:)),'k')
    hold on
    plot(mean(response(bgClassX==2,:)),'b')
    plot(mean(response(bgClassX==1,:)),'r')
    
% makeAxisStruct(gca,strtrim(['psthMeans' experimentDate cellNum]))

frequencyStrip=seqMean; %freqs contained in response

fs=10000; %sampling frequency
L=length(frequencyStrip);
NFFT = 100000;
X = fftshift(fft(frequencyStrip,NFFT));
Pxx=X.*conj(X)/(NFFT*NFFT); 
f = fs*(-NFFT/2:NFFT/2-1)/NFFT; %Frequency Vector
figure(10)
plot(f,abs(X)/(L),'r');
title('Magnitude of FFT');
xlabel('Frequency (Hz)')
ylabel('Magnitude |X(f)|');
xlim([0 100])
% makeAxisStruct(gca,strtrim(['powerSpectrum' experimentDate cellNum]))

% X = fft(frequencyStrip,NFFT);
% X = X(1:NFFT/2+1);%Throw the samples after NFFT/2 for single sided plot
% Pxx=X.*conj(X)/(NFFT*NFFT);
% f = fs*(0:NFFT/2)/NFFT; %Frequency Vector
% figure(12)
% plot(f,10*log10(Pxx),'r');
% title('Single Sided - Power Spectral Density');
% xlabel('Frequency (Hz)')
% ylabel('Power Spectral Density- P_{xx} dB/Hz');
% xlim([0 100])
% 
% fPower = f';
% logY=(10*log10(Pxx))';


%% begin LN - Grab Frames

noiseVars = struct();
noiseVars.type = gaussianOrBinary;
noiseVars.contrast = 0.3333;
timings = [250,10000,250]; % AUTOMATE THIS LATER
frames = manookinlab.ovation.getFrameTimesFromMonitor(frameTimings(usableIndex,:), 10000, binRate);
frameValues = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,2,double(string(seedList(usableIndex)))',frameDwell(usableIndex));

%% Generate Filter & NL (incl. modeled curves)
% 

colorSet = [1 0 0;0 0 1;0 0 0];
lfilter=[]; 
nonlinearityBins=100;
   stimDuration = preTime + (1 : stimTime);     
for f = 1:3      
    stimulus = frameValues(bgClassX==f, :);
    responses = response(bgClassX==f, :);

    lfilter(f,:) = getLinearFilter(stimulus, responses, ... %
        'analysisType', 'revcorr', ...
        'fourierCorrection',false, ...
        'binRate', binRate, ...
        'filterTime', 0.5, ...
        'frameRate', frameRate);
    lfilter(f,:) = lfilter(f,:) / norm(lfilter(f,:));
end
allFilters = lfilter;
lfilter = mean(lfilter,1); %without mean filter, NLs sometimes look different
figure(5)
plot(lfilter(1:500),'m','LineWidth',2)
hold on
plot(allFilters(1,1:500),'--r')
plot(allFilters(2,1:500),'--b')
plot(allFilters(3,1:500),'--k')
legend('AverageFilter','SequentialFilter','RandomFilter', 'StationaryFilter')

title('sta')
xlabel('ms')
ylabel('contrast weight')
for g = 1:3
    
    % nl starts here
    Stim = frameValues(bgClassX==g,:);
    Resp = response(bgClassX==g,:);
    lf = lfilter(1 : size(Stim,2));
    Pred = zeros(size(Resp));
        for p = 1 : size(Stim,1)
            revcor = ifft( fft(Stim(p,:)) .* fft(lf) );
            Pred(p,:) = real(revcor);
            
            Pred(p,:) = Pred(p,:)./std(Pred(p,:));
        end
    [xBin(g,:),yBin(g,:)] = binNonlinearity(Pred(:,stimDuration(501:end)),Resp(:,stimDuration(501:end)),nonlinearityBins);
%  [xBin(f,:),yBin(f,:)] = binNonlinearity(Pred(:,stimDuration(:)),Resp(:,stimDuration(:)),nonlinearityBins);

    figure(1)
nlParams(g,:) = fitNonlinearityParams(xBin(g,:), yBin(g,:));
% nlX(f,:) = linspace(-max(abs(xBin(f,:))),max(abs(xBin(f,:))),10000);
% plot(nlX(f,:),outputNonlinearity(nlParams(f,:),nlX(f,:)))    
plot(xBin(g,:),outputNonlinearity(nlParams(g,:),xBin(g,:)))    
hold on

figure(2)
plot(xBin(g,:),yBin(g,:),'Color',colorSet(g,:))
hold on
    legend('Sequential','Random','Static')
    ylabel('Spike Rate (Hz)')
end

% makeAxisStruct(gca,strtrim(['rawishNLs' experimentDate cellNum]))
%% Verify Model -- go back to sorting code to get repeated seeds

noiseVars = struct();
noiseVars.type = gaussianOrBinary;
noiseVars.contrast = 0.3333;
timings = [250,10000,250]; % AUTOMATE THIS LATERaa
frames = manookinlab.ovation.getFrameTimesFromMonitor(frameTimings(usableIndex,:), 10000, binRate);
frameValues = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),1000,frames,2,double(string(seedList(usableIndex)))',frameDwell(usableIndex));
bgClassRepSeeds = bgClass(usableIndex);
responseRepSeed = psthMatrix(usableIndex,:);
%% Verify Model -- generate model estimate
avgFilter=mean(lfilter(:,1:size(frameValues,2)));

rSeedGenerator = ifft(fft(avgFilter) .* fft(frameValues(1,:)));
% rSeedGenerator = real(rSeedGenerator);
rSeedGenerator = rSeedGenerator/std(rSeedGenerator);
rSeedY = outputNonlinearity(nlParams(1,:),rSeedGenerator);
plot(responseRepSeed(1,:));
hold on 
plot(rSeedY);
title('off smooth - sequential')
legend('data','model')
xlabel('input')
ylabel('Spike Rate (Hz)')
