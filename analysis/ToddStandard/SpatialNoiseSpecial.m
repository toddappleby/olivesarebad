%% Spatial  Noise

exportDirectory = 'C:\Users\reals\Documents\PhD 2021\ClarinetExports\';
experimentDate = '2023_0808';
cellNum = 'Ac4';
cd(strcat(exportDirectory,experimentDate))

cd(strcat(exportDirectory,experimentDate)) 
load(strcat(experimentDate(1:4),experimentDate(6:9),cellNum,'_FT.mat'))
frameTs = epochs;
load(strcat(experimentDate(1:4),experimentDate(6:9),cellNum,'.mat'))
cellData = epochs;

%%

% splitFactors = ["frameDwell","noiseClass"];

% splitFactors = ["frameRate","chromaticClass","contrast"];


splitFactors = ["stixelSize","chromaticClass","contrast"];
% epochGroup = 'Split1';
epochGroup = 'Control';
% 

desiredSTD = 5;

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end


count = 0;
clear epochStorage;
count2 =0;


for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
    egLabel = epochs(i).meta.epochGroupLabel;
    if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
    end
   
%                   if strcmp(displayName,'Spatial Noise') && strcmp(egLabel,'Control')  
                  if strcmp(displayName,'Jittered Noise') && strcmp(egLabel,epochGroup) 
%              if strcmp(displayName,'Fast Noise') && strcmp(egLabel,epochGroup)  
                for s=1:length(splitCell)  
                  if strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'double')
                      stringedEntry = convertCharsToStrings(num2str(getfield(epochs(i).meta,splitFactors(s))));
                  elseif strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'char')
                      stringedEntry = convertCharsToStrings(getfield(epochs(i).meta,splitFactors(s)));
                  end
               
                    splitCell{2,s}=[splitCell{2,s} stringedEntry];
                end

            count = count + 1;
%             intensity(count) = epochs(i).meta.intensity;
            epochStorage(count,:) = epochs(i).epoch;
            
            frameTimes(count,:) = frameTs(i).epoch;
            numXChecks(count) = epochs(i).meta.numXChecks;
            numYChecks(count) = epochs(i).meta.numYChecks;
            numXStixels(count) = epochs(i).meta.numXStixels;
            numYStixels(count) = epochs(i).meta.numYStixels;
            numFrames = epochs(i).meta.numFrames;
            stepsPerStixel = epochs(i).meta.stepsPerStixel;
            seeds(count)= epochs(i).meta.seed;
            preTime = epochs(i).meta.preTime;
            stimTime = epochs(i).meta.stimTime;
            tailTime = epochs(i).meta.tailTime;
              end

              
end
       
count = 0;
clear monitorStorage
clear epochStartTime



frames = manookinlab.ovation.getFrameTimesFromMonitor(frameTimes, 10000, 1000);
        
      
 


stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),(size(epochStorage,2))/10);

%string thing
% for m = 1:size(splitCell,2)-1
%     if strcmp(class(splitCell{2,m}),'char')
%         uniqueStrings = unique(splitCell{2,m});
%         replacer = zeros(size(splitCell{2,m}));
%         
%         for n = 1:length(uniqueStrings)
%            strInd = find(splitCell{2,m}==uniqueStrings(n));
%            replacer(strInd) = n;
%         end
% %         holdList(:,m) = splitCell{2,m};
%         splitCell{2,m}=replacer;
%     end
% end

allSets=strings(size(splitCell{2,1},2),size(splitCell,2)-1);

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
%         psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,10);
        psthMatrix(k,:) = psth(spikeMatrix(k,:),10,sampleRate,10);
%         psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,1);
          end
end

%60 f/sec * 21sec 
% numFrames = 1260;
binnedSpike=[];
for s = 1:size(spikeMatrix,1)
binnedSpikes(s,:) = binSpikeCount(spikeMatrix(s,:),1000,10000)';
end

params = struct();
params.saveGraph =0;
params.stimName = 'bars';
timings = [preTime stimOrig tailTime];
% [allM,frameVals] = doSpatialMap(psthMatrix,numXChecks(1),numYChecks(1),numFrames,seeds,timings);
[allM,allData,fvYellow,fvBlue,frameValsCombined,frameTs,t1,t2] = doJitteredMap(binnedSpikes(indexHolder{2,1},:),max(numXChecks),max(numYChecks),max(numYStixels),max(numXStixels),indexHolder{1,1}(2),numFrames,seeds(indexHolder{2,1}),timings);
% [allM,allData,frameVals,frameValsCombined,frameTs,t1,t2] = doJitteredMap(binnedSpikes(indexHolder{2,2},:),min(numXChecks),min(numYChecks),min(numYStixels),min(numXStixels),indexHolder{1,2}(2),numFrames,seeds(indexHolder{2,2}),timings);
% [allMBlue,allMYellow,allData,fvBlue,frameValsCombined,frameTs,t1,t2,fvYellow] = doFastMap(binnedSpikes,numXChecks(1),numYChecks(1),numYStixels(1),numXStixels(1),numFrames,seeds,timings);
%% Interspike Interval

allTimes = [];
isiTmp=[];

for t = 1:size(spikeMatrix,1)
    
    spikeTimes = find(spikeMatrix(t,:)>0);
    isiTmp = diff(spikeTimes/10)'; %samplerate / binrate for milleseconds
    allTimes = [allTimes; isiTmp];
    
     
end

[y,x]= getSpikeAutocorrelation(allTimes); 
 
figure
plot(y)
%% plot test 1

for c = 1:size(allM,3)
    imagesc(allM(:,:,c))
    
    pause
end

% for c = 1:size(sta,3)
%     imagesc(sta(:,:,c))
%     pause 
% end
% figure(5) 
% for p = 1:6        
% subplot(2,3,p)
% imagesc(sta(:,:,p+1))
% title(p+1)
% end


%% RGB Matrix

% RGBMatrix = zeros(size(fvYellow,1),size(fvYellow,2),size(fvYellow,3),3);
% 
% RGBMatrix(:,:,:,1)=fvYellow;
% RGBMatrix(:,:,:,2)=fvYellow;
% RGBMatrix(:,:,:,3)=fvBlue;



RGBMatrix(:,:,1) = fvYellow(:,:,1);
RGBMatrix(:,:,2) = fvYellow(:,:,1);
RGBMatrix(:,:,3) = fvBlue(:,:,1);

imagesc(RGBMatrix)
figure
imagesc(fvYellow(:,:,1))
title('yellow')
colormap('gray')
figure
imagesc(fvBlue(:,:,1))
title('blue')
colormap('gray')
%% try single epoch with actual frames


for c = 1:size(binnedSpikes,1)  %number of epochs

    binsPerFrame = 1 ; 
    fTimes = frames{c};

    spikeArray = binnedSpikes(c,:);

    firstFrame = find(fTimes<timings(1),1,'last')+1;
%     firstFrame=16;
    
if c==1
    lastFrame = find(fTimes>(timings(2)+timings(1)),1,'first');
end
%     totalFrames = (lastFrame - firstFrame) +1;
    %fTimes(firstFrame:lastFrame);

    [frameValues,fvBlue,~,fvAll] = getFastNoiseFrames(numXStixels(1),numYStixels(1), numXChecks(1), numYChecks(1), 'BY', lastFrame, 2, seeds(c));
%     frameValues = getJitteredNoiseFrames(numXStixels(1),numYStixels(1), numXChecks(1), numYChecks(1), lastFrame, 2, seeds(c));


    bdata = zeros(1,lastFrame);

     for n = 1 : lastFrame %binPerFrame = 1
            timePerFrame = fTimes(firstFrame+(n-1)) : fTimes(firstFrame+n)-1; %-1 because span until just before next frame
            bdata(n) = mean(spikeArray(1,timePerFrame));
     end

    bdata = bdata - median(bdata(60*binsPerFrame+1:end));
        bdata(1:60*binsPerFrame) = 0;
        
% size(revCorMatrix)
        if c == 1
            
    revCorMatrix = zeros(size(frameValues,1),size(frameValues,2),lastFrame);
        end
    for m = 1:size(frameValues,1)
        for n = 1:size(frameValues,2)   
            revCorMatrix(m,n,:) = squeeze(revCorMatrix(m,n,:)) + fft(bdata(:)) ...
                    .* conj(fft(squeeze(frameValues(m,n,:))));
        end
    end

end

sta = zeros(numYChecks(1),numXChecks(1),1276);
    for k = 1 : numYChecks(1)
        for m = 1 : numXChecks(1)
            tmp = ifft(squeeze(revCorMatrix(k,m,:)));
%              plot(tmp(1:30))  
%              pause
            sta(k,m,:) = tmp(1 : 1276);
        end
    end

%% LN NL LN NL LN 

% stim = frameValues(17,30,:);
lf = sta(22,26,:);
lf = squeeze(lf)/norm(squeeze(lf));
for z = 1:size(binnedSpikes,1)

    fTimes = frames{z};
    spikeArray = binnedSpikes(z,:);
    firstFrame = find(fTimes<timings(1),1,'last')+1;
    lastFrame = find(fTimes>(timings(2)+timings(1)),1,'first');
    
    frameValues = getJitteredNoiseFrames(numXStixels(1),numYStixels(1), numXChecks(1), numYChecks(1), lastFrame, 2, seeds(z));

    
     for n = 1 : lastFrame
            timePerFrame = fTimes(firstFrame+(n-1)) : fTimes(firstFrame+n)-1; %-1 because span until just before next frame
            bdata(n) = mean(spikeArray(1,timePerFrame));
     end
     
     bigData(z,:) = bdata;
     bigStim(z,:) = frameValues(22,26,:);

end

   Pred = zeros(size(bigData));
     for p = 1 : size(bigStim,1)
            revcor = ifft( fft(bigStim(p,:)) .* fft(squeeze(lf)'));
            Pred(p,:) = real(revcor);
            
            Pred(p,:) = Pred(p,:)./std(Pred(p,:));
%             Pred(p,:) = Pred(p,:);
     end
Resp=bigData;
Resp=Resp*1000;
[xBin,yBin] = binNonlinearity(Pred(:,17:end),Resp(:,17:end),100);


figure(9)
hold on
plot(xBin,yBin)
%% LN NL LN NL LN 

% stim = frameValues(17,30,:);

gridY = 31; %bottom on
gridX = 38;
% 
% gridY = 16;
% gridX = 32; %top On
% 
% gridY = 23;
% gridX = 37;  %off center

lf = sta(gridY,gridX,:);
lf = squeeze(lf)/norm(squeeze(lf));
figure(2);hold on;plot(lf(1:30),'Color',[0.9290 0.6940 0.1250],'LineWidth',2) %[0.6350 0.0780 0.1840]
% figure(2);hold on;plot(lf(1:30),'Color',[0.6350 0.0780 0.1840],'LineWidth',2) %

for z = 1:size(binnedSpikes,1)

    fTimes = frames{z};
    spikeArray = binnedSpikes(z,:);
    firstFrame = find(fTimes<timings(1),1,'last')+1;
    lastFrame = find(fTimes>(timings(2)+timings(1)),1,'first');
    
    frameValues = getJitteredNoiseFrames(numXStixels(1),numYStixels(1), numXChecks(1), numYChecks(1), lastFrame, 2, seeds(z));

    
     for n = 1 : lastFrame
            timePerFrame = fTimes(firstFrame+(n-1)) : fTimes(firstFrame+n)-1; %-1 because span until just before next frame
            bdata(n) = mean(spikeArray(1,timePerFrame));
     end
     
     bigData(z,:) = bdata;
     bigStim(z,:) = frameValues(gridY,gridX,:);

end

   Pred = zeros(size(bigData));
     for p = 1 : size(bigStim,1)
            revcor = ifft( fft(bigStim(p,:)) .* fft(squeeze(lf)'));
            Pred(p,:) = real(revcor);
            
            Pred(p,:) = Pred(p,:)./std(Pred(p,:));
%             Pred(p,:) = Pred(p,:);
     end
Resp=bigData;
Resp=Resp*1000;
[xBin,yBin] = binNonlinearity(Pred(:,17:end),Resp(:,17:end),100);

figure(7)
hold on
nlParams = fitNonlinearityParams(xBin, yBin);   
plot(xBin,outputNonlinearity(nlParams,xBin))    


figure(9)
hold on
plot(xBin,yBin)


%% STACK FRAMES 

m = 49;
numFrames = 1815;
stackedFrames = zeros(1,numFrames*m);
for z = 1:m
    
    frameHold = frames{z};
    frameHold = frameHold(find(frameHold>250,1):end);
    frameHold((numFrames+1):end)= [];
    
    stackedFrames((1+(numFrames*(z-1))):(numFrames+(numFrames*(z-1)))) = frameHold+(21750*(z-1));
    
end
%% STACK SPIKE TIMES
spTimes = [];
m = 93;
for s = 1:m
spTimeHolder = find(spikeMatrix(s,:)==1);
spTimeHolder = spTimeHolder/10;
% if ~isempty(spTimes)
spTimes = [spTimes (spTimeHolder+(21750*(s-1)))];
% else
%     spTimes = [spTimes spTimeHolder];
% end
end
%% save important variables for python 
save('importStimforPython.mat','allM','frameValsCombined','stackedFrames','spTimes','allData');


%% reformat frame vals for python

% for g = 1:size(frameVals,3)
%     stimFrames(g,:,:)=frameVals(:,:,g);
% end
binRate=1000;

time = linspace(0,21,(21000*size(psthMatrix,1)));
psthFull = psthMatrix(:,251:21250);
psthFull = psthFull(:)';
% spikesFull = spikeMatrix(:)';


% spikes = [];
% for d = 1:size(spikeMatrix,1)
%     spikes(d,:) = binSpikeCount(spikeMatrix(d,2501:212500), binRate, sampleRate)';
%   
%     
% end
% spikesFull = spikes';
% spikesFull= spikesFull(:);
% spTimes = find(spikesFull>0)';
%  spTimes = find(allData > 0);
spTimes = find(spikeMatrix(1,:)==1);
spikesInOrder = allData;

spTimes2 = find(spikesInOrder>0);
% spTimes = spTimes/10;%ms

psthAllData=psth(allData,10,sampleRate,1);
% frameTs=frameTs/10; %ms
spTimes=spTimes/10;
frameTs = frames{1};
% frameTs = frameTs/1000;
% frameTs=frameTs(find(frameTs<21250),1)
frameTs=frameTs(find(frameTs>250,1):end)
%THEN, take out first 16 frame values
%OR, do 1275 and push spikes times back 250 ms.... 
% frameTs(1276:end)=[];
%save('importStimforPython.mat','allM','frameValsCombined','frameTs','time','spTimes','allData','psthAllData');
save('importStimforPython.mat','allM','frameValsCombined','frameTs','time','spTimes','spTimes2','allData','psthAllData');

%stim formatting 


%% comparing spike times to frame times

timeAxis = zeros(1,21500);
for k = 1:length(spTimes)
timeAxis(round(spTimes(k))) = 1;
end

timeAxis2= zeros(1,21500);
for l = 1:length(frameTs)
    timeAxis2(round(frameTs(l))) = 1;
end

%% frame vals?? timing
% should be 1260 framse whilte framerate and stimTime remain unchanged stim
% parameters ]]
howmany = round(21000/16.6667);
% frameFormatter = zeros(size(spikes,1),howmany); 
frameFormatter = [];
frameMS=17;
totalFrames = 1260 * size(spikes,1);
for f = 1:totalFrames
        for p = 1:frameMS
frameFormatter(:,:,f)=frameValsCombined(:,:,p+(frameMS*(f-1)));
        end
end


%% plot test python  


for q = 1:size(frameVals,3)
    imagesc(frameVals(:,:,q))
    pause
end


%% plot test 1

for c = 1:size(allM,3)
    imagesc(allM(:,:,c)) 
    pause
end
%%

for z = 4:2:18
   subplot(1,11,z-3)
   imagesc(allM(:,:,z))

end


%% plot test 2
for p = 5:16
subplot(4,3,p-4)
imagesc(allM(:,:,p))
title(p+1)
end
%% test
 tester = allM(24,24,:);
tester2 = ones(1,size(t1,1)).*tester(:);
normStimSet = t1(:,1:30)-tester2';
% normStimSet = t1(:,1:30);


covarM = (1/29)*normStimSet' * normStimSet;
[U S V] = svd(covarM);

%% rank thing
rank = 1;

test = U(:,1:rank)*S(1:rank,1:rank)*V(:,1:rank)';
figure
plot(test')
% figure(2)
% colormap(gray)
% imagesc(reshape(test(30,:),10,10))
%% 
spikeStims = t1(:,1:30);
% repmat(mn,1,n)
[m,n]=size(spikeStims);   
mn=mean(spikeStims, 1); %sta (basically)
spikeStims=spikeStims-mn; % centering stims
Cx=(1/(n-1))*spikeStims'*spikeStims;  % Create covariance matrix
[V,D]=eig(Cx); %run decomposition to find eigenvalues & eigenvectors
lambda=diag(D);  % get eigen values from matrix diagonals

[sortthing,m_arrange]=sort(-1*lambda);  %eigenvalues come out backwards 
V=V(:,m_arrange);  

Y=V'*Cx; % project covariance matrix along eigenvector axis

sortthing = -sortthing;
plot(sortthing/max(sortthing))
title('Distribution of Mode Energy from STC Analysis')
xlabel('Eigenvalue Number')
ylabel('Normalized Energy')


