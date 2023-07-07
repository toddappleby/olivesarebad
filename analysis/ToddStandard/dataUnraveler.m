%% load FT first, give different name --- #1
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2019_1107') %2019 1107 bc7 sbc
cellDate = '2021_0122';
cellNum = 'Ac1';
cellName = strcat(cellDate,cellNum);

loadDate = char(cellDate);
loadDate = strcat(loadDate(1:4),loadDate(6:9));

fileName = strcat(loadDate,cellNum,'.mat');             

% directory = string("/Users/toddappleby/Documents/Data/Clarinet Exports/" + cellDate);
directory = string(strcat(getDirectory,cellDate)); %WINDOWS <<

cd(directory)

load(fileName)
% frameTimings = epochs;

% get protocols 
for z = 1:size(epochs,2)
   list(z) = string(epochs(z).meta.displayName);
end
uniqueProtocols = unique(list')
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
% saveStuff.flag = saveFlag;
% saveStuff.cellType = cellType;
% saveStuff.expNum = expNum;
% saveStuff.OFFParasolC = OFFParasolC;
% saveStuff.ONParasolC = ONParasolC;

logi = contrast == .05;
[cSensitivity, ~, ~] = frequencyModulation(mean(psthMatrix(logi,5001:end ),1), 10000, 4, 'avg', 1, []);
%% chirp regurgitation
dataType = 1;
clear epochStorage
saveFlag = 0;
count =0;
maxFreq = 25; %at least 1 has 22 randomly ... maybe others?

controlOrAlt = 'Control';
% controlOrAlt = 'AltCenter';
drawStruct = struct(); 

desiredSTD=4; 

if dataType == 1
for i = 1:length(epochs)
   displayName = epochs(i).meta.displayName;

   recordingTechnique = epochs(i).meta.recordingTechnique;
  egLabel = epochs(i).meta.epochGroupLabel;
%    
   
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
   end
   
    if isfield(epochs(i).meta,'led')
       ledTHING = epochs(i).meta.led;
   end



   if strcmp(displayName,'Chirp Stimulus') && strcmp(oAnalysis,'extracellular') && epochs(i).meta.frequencyMax == maxFreq && strcmp(egLabel,controlOrAlt)
count = count + 1;
epochStorage(count,:) = epochs(i).epoch;


% contrast(count) = epochs(i).meta.contrast;

preTime = epochs(i).meta.preTime;

stimTime = epochs(i).meta.stimTime;

tailTime = epochs(i).meta.tailTime;

drawStruct(count).frequencyMax = epochs(i).meta.frequencyMax;
drawStruct(count).frequencyMin = epochs(i).meta.frequencyMin;
drawStruct(count).frequencyTime = epochs(i).meta.frequencyTime;
drawStruct(count).interTime = epochs(i).meta.interTime;
drawStruct(count).frequencyContrast = epochs(i).meta.frequencyContrast;
drawStruct(count).contrastMax = epochs(i).meta.contrastMax;
drawStruct(count).contrastMin = epochs(i).meta.contrastMin;
drawStruct(count).contrastTime = epochs(i).meta.contrastTime;
drawStruct(count).contrastFrequency = epochs(i).meta.contrastFrequency;
drawStruct(count).preTime = epochs(i).meta.preTime;
drawStruct(count).tailTime = epochs(i).meta.tailTime;
drawStruct(count).stepTime = epochs(i).meta.stepTime;
drawStruct(count).sampleRate = epochs(i).meta.sampleRate;
drawStruct(count).stepContrast = epochs(i).meta.stepContrast;

   end

end
else
  for i = 1:length(epochs)
   displayName = epochs(i).meta.displayName;

   recordingTechnique = epochs(i).meta.recordingTechnique;

   egLabel = epochs(i).meta.epochGroupLabel;

   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
   end

   if strcmp(displayName,'Chirp Stimulus') && strcmp(oAnalysis,'analog') 
count = count + 1;
epochStorage(count,:) = epochs(i).epoch;

% contrast(count) = epochs(i).meta.contrast;
% 
% preTime = epochs(i).meta.preTime;
% 
% stimTime = epochs(i).meta.stimTime;
% 
% tailTime = epochs(i).meta.tailTime;

   

 

        end

    end  
end

for z = 1:size(epochStorage,1)
figure(1)
plot(epochStorage(z,:))
set(gca,'visible','off')

hold on

end
clear spikeMatrix
clear psthMatrix
sampleRate=10000;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000
if dataType==1
for k = 1:size(epochStorage,1)
    
[spikes, finalSTD, finalDiscard] = convertSpikesAdree(epochStorage(k,:), stimTime, ...
              desiredSTD);
          if isempty(spikes)
              disp('deleted epoch')
          else
        spikeMatrix(k,:) = spikes;
        psthMatrix(k,:) = psth(spikeMatrix(k,:),10,sampleRate,1);
          end
end
else
    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
    end
end

figure(3)
xvalsChirp = linspace(1,23000,230000);
plot(xvalsChirp,mean(psthMatrix,1)')

figure(2)
chirpStim = drawChirp(drawStruct);
plot(mean(epochStorage,1)/max(mean(epochStorage,1))+3)
hold on
plot(chirpStim)

figure(4)
plot(mean(psthMatrix,1)/max(mean(psthMatrix,1))+2)
hold on
plot(chirpStim)

figure(7)
plot(mean(spikeMatrix,1)/max(mean(spikeMatrix,1))+1.1)
hold on
plot(chirpStim)

freqDuration = drawStruct(1).frequencyTime * 10; %freq sweep in time points ! 
%*** Frequency sweep stuff.  There might be more to do here.  
% chirpFFT(psthMatrix,drawStruct);


%**** Contrat sweep -- problematic if length of sweep changes between
%epochs.  I'll have to switch to the normal sorting stuff eventually
% meanRate = mean(mean(psthMatrix(:,140000:145000),2),1);
% contrastPeaks = findpeaks(chirpStim(1,145000:225000));
% 
% contrastPeaks = 2*(contrastPeaks-.5); % contrast relative to bg
% contrastResponse = smoothdata(mean(psthMatrix(:,145000:225000),1)); %only looking for peak effect,so smoothing fine--actually might make peak rate more realistic becuz binning makes peak v high
% contrastResponsePeaks = findpeaks(contrastResponse,'MinPeakDistance',2500); %find data peaks from smoothed contrast resp data
% contrastResponsePeaks = contrastResponsePeaks - meanRate; %subtract mean rate for low contrast effect
% figure(8)
% if length(contrastPeaks) == length(contrastResponsePeaks)
% plot(contrastPeaks,contrastResponsePeaks); %probably better way to ensure matching peaks, but above lines should make them match anyway
% title('Contrast Sweep: Peak Response at each Cycle')
% xlabel('Contrast (relative intensity)')
% ylabel('Peak Spike Rate (Hz)')
% else
%     peakDiffy = abs(length(contrastResponsePeaks) - length(contrastPeaks));
%         if length(contrastResponsePeaks) > length(contrastPeaks)
%             
%             plot(contrastPeaks,contrastResponsePeaks((peakDiffy+1):end)); %remove initial peaks (usually these are the problem...)
%             disp('WARNING: MISMATCHED PEAKS, SUGGEST MANUAL CHECK OF PEAKS')
%         elseif length(contrastResponsePeaks) < length(contrastPeaks)
%             contrastResponsePeaks = [zeros(1,peakDiffy) contrastResponsePeaks];
%             plot(contrastPeaks,contrastResponsePeaks)
%         end
% end
%***

if saveFlag == 1 
     maxFreq = drawStruct(size(psthMatrix,1)).frequencyMax;
     meanChirp = mean(psthMatrix,1)';
     xvalsChirpSave = xvalsChirp';
%      saveName = strcat('chirp_',cellName);
contrastPeaks = contrastPeaks';
contrastResponsePeaks = contrastResponsePeaks';
     save('Chirp.mat','meanChirp','xvalsChirpSave','maxFreq')
     save('chirpContrast.mat','contrastPeaks','contrastResponsePeaks')
end
% chirp tests

freqStart = (drawStruct(1).preTime + drawStruct(1).stepTime*2 + drawStruct(1).interTime*2)*10;%pts
freqEnd = freqStart + drawStruct(1).frequencyTime * 10; %pts

freqDeriv = diff(chirpStim(freqStart:freqEnd)); %velocity

freqSweep = chirpStim(freqStart:freqEnd);

chirpSpikes = mean(spikeMatrix(:,freqStart:freqEnd),1); 

chirpSpikesOff = chirpSpikes(freqDeriv <0);
chirpSpikesOn = chirpSpikes(freqDeriv>0);

offLogical = freqDeriv < 0;
offVector = find(offLogical ==1);
onLogical = freqDeriv>0;
onVector = find(onLogical ==1);
g=1;
h=1;
offCycleSpikes = [];
onCycleSpikes = [];
escapeVectorOff=[];
escapeVectorOn=[];
for t = 2:length(offVector)
    if offVector(t) - offVector(t-1) == 1
        escapeVectorOff(g,t-1) = offVector(t-1); % works because first cycle always biggest
    else
        g=g+1;
    end
    

end
    
for t = 2:length(onVector)
    if onVector(t) - onVector(t-1) == 1
        escapeVectorOn(h,t-1) = onVector(t-1); % works because first cycle always biggest
    else
        h=h+1;
    end
end

for b = 1 : g
   
     offOut = escapeVectorOff(b,escapeVectorOff(b,:)>0);
     cutOut1 = floor(.05*length(offOut));
     onOut = escapeVectorOn(b,escapeVectorOn(b,:)>0);
     cutOut2 = floor(.05*length(onOut));
     
     offCycleSpikes(1,b)=sum(chirpSpikes(1,offOut(cutOut1:end)),2);
     onCycleSpikes(1,b)=sum(chirpSpikes(1,onOut(cutOut2:end)),2);
end

% for c = 1:h
freqXVals = linspace(drawStruct(1).frequencyMin,drawStruct(1).frequencyMax,b);

figure(9)
plot(freqXVals,offCycleSpikes,'r')
hold on

plot(freqXVals,onCycleSpikes,'b')


legend('off response','on response')
xlabel('frequency sweep (Hz)')
ylabel('mean spike count')

%hill

coef = [1 -2];
% fitcoef = nlinfitsome([false false], freq, amp ./ max(amp), @hill, coef);
% fit = hill(fitcoef, freq);
freqMax = drawStruct.frequencyMax;
tester = linspace(0,freqMax,1000);

maxOff = max(offCycleSpikes);
maxOn = max(onCycleSpikes);
if maxOn > maxOff 
    max1 = maxOn;
    max2 = 1;
else
    max2 = maxOff;
    max1 = 1;
    
end


fitcoef = nlinfitsome([false false], freqXVals,offCycleSpikes./max(offCycleSpikes),@hill,coef);
fit = hill(fitcoef,tester);
figure(51)
plot(tester,fit/max1,'r')
hold on
fitcoef = nlinfitsome([false false], freqXVals,onCycleSpikes./max(onCycleSpikes),@hill,coef);
fit = hill(fitcoef,tester);
plot(tester,fit/max2,'b')



%% MTF spots and annuli
dataType=1;
% currWanted='excitation';
currWanted='inhibition';
splitFactors = ["contrast","stimulusClass","temporalFrequency","temporalClass","chromaticClass","radius"];
% splitFactors = ["temporalFrequency","searchAxis","barSize","position"];
runSplitter = ["radius"]; 

controlOrAlt =  "Control";
% controlOrAlt =  "NDF 4";
% controlOrAlt =  "Wholecell_exc";
% controlOrAlt =  "AltCenter";

splitCell = cell(2,length(splitFactors));

desiredSTD =6;

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage 

stringComparer = [];
stringComparer = string(stringComparer);
if dataType ==1
for i = 1:length(epochs)
    
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
%    recording = epochs(i).meta.recordingTechnique;
  
   
  if isfield(epochs(i).meta,'onlineAnalysis')
      oAnalysis = epochs(i).meta.onlineAnalysis;
  end
  

   

%    if strcmp(displayName,' Fspot') 
   if strcmp(displayName,'S MT Fspot') && strcmp(egLabel,controlOrAlt)
    temporalClass = epochs(i).meta.temporalClass;
%     if strcmp(temporalClass,'sinewave')
       if strcmp(temporalClass,'pulse') 
        count = count + 1;
        for s = 1:length(splitFactors)
            if ischar(getfield(epochs(i).meta,splitFactors(s)))
                splitCell{2,s} = [splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];    
            else
                splitCell{2,s}= [splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
            end
        end
        
      
        
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
end

else
    for i = 1:length(epochs)
    
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
%    recording = epochs(i).meta.recordingTechnique;
  
   
  if isfield(epochs(i).meta,'onlineAnalysis')
      oAnalysis = epochs(i).meta.onlineAnalysis;
  end
  
currType=[];
          if strcmp(egLabel,'Whole cell_exc') || strcmp(egLabel, 'Wholecell_exc')
           currType = 'excitation';
       elseif strcmp(egLabel,'Whole cell_inh') || strcmp(egLabel, 'Wholecell_inh')
           currType = 'inhibition';
       end

%    if strcmp(displayName,' Fspot') 
   if strcmp(displayName,'S MT Fspot') && strcmp(currWanted,currType)
    temporalClass = epochs(i).meta.temporalClass;
    if strcmp(temporalClass,'sinewave')
%        if strcmp(temporalClass,'pulse')
        count = count + 1;
        for s = 1:length(splitFactors)
            if ischar(getfield(epochs(i).meta,splitFactors(s)))
                splitCell{2,s} = [splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];    
            else
                splitCell{2,s}= [splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
            end
        end
        
      
        
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
 end
end
stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));


% replacer = zeros(size(splitCell{2,1));



for m = 1:size(splitCell,2)-1
    if strcmp(class(splitCell{2,m}),'string')
        uniqueStrings = unique(splitCell{2,m});
        replacer = zeros(size(splitCell{2,m}));
        
        for n = 1:length(uniqueStrings)
           strInd = find(splitCell{2,m}==uniqueStrings(n));
           replacer(strInd) = n;
        end
%         holdList(:,m) = splitCell{2,m};
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

if dataType == 1
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

else


    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
    end
end

timings = [preTime stimOrig tailTime];
params = struct();
params.saveGraph = 0;
params.stimName = 'mTF';
params.tfreq = unique(tfreq);
params.dataType = dataType;
params.cellName= cellName;
% saveGraph = 0;
% stimName = 'Grating';
% freqAnalysis(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
simpleSpots(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% 
plot(excitation(13,:))
hold on
plot(inhibition(13,:)*1.38)

scaledInhibition = inhibition.*1.38;

diffCurrents = scaledInhibition-excitation;

meanOffCurrents = mean(diffCurrents(:,5001:10000),2);
meanOffCurrents(1)=0;
plot(unique([splitCell{2,6}]'),meanOffCurrents)
hold on
meanOnCurrents = mean(diffCurrents(:,10001:15000),2);
meanOnCurrents(1)=0;
plot(unique([splitCell{2,6}]'),meanOnCurrents)
%% Bar centering 
desiredSTD = 5;
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
   epochLabel = epochs(i).meta.epochGroupLabel;
if strcmp(displayName,'Bar Centering') && strcmp(epochLabel,'Control')
   
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
%% Simple Response (Color spot and OMS)
dataType = 1; 
desiredSTD = 5;
% % 
splitFactors = ["stimulusClass"];
% splitFactors = ["radius","splitContrasts","contrast","spaceConstant"];
% splitFactors = ["stimulusIndex"];
% splitFactors = ["outerRadius","chromaticClass","contrast"];
% splitFactors = ["chromaticClass"];
% splitFactors = ["led"];
% splitFactors = ["backgroundIntensity","spotIntensity"];
% %
% protocolID = "Doves Movie";
protocolID = "Object Motion Grating";
% protocolID = "Object Motion Dots";
% protocolID = "Single Spot";
% protocolID = "Chromatic Spot";
% protocolID = "Led Pulse";
splitCell = cell(2,length(splitFactors));


for g = 1:size(splitCell,2)
    splitCell{1,g} = splitFactors(g);
end

controlOrAlt = "Control";
% controlOrAlt = "Whole cell_exc";
% controlOrAlt = "rf annuli";
% controlOrAlt = "AltCenter";
% controlOrAlt = "AltCenter_exc"

  counter = 0;
  epochLength =0;
for a = 1:length(epochs)
   displayName = epochs(a).meta.displayName;
   falseLength = length(epochs(a).epoch);
   if strcmp(displayName,protocolID) 
%    if strcmp(displayName,'Chromatic Spot')
%        if strcmp(displayName,'Led Pulse')
       counter = counter +1;
       
           if falseLength > epochLength
               finale = falseLength;
           end
           
       epochLength = length(epochs(a).epoch);
   end
end

epochStorage = zeros(counter,finale);
count = 0;
% clear epochStorage;
if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   egLabel = epochs(i).meta.epochGroupLabel;
%    if strcmp(displayName,'Chromatic Spot') && ~strcmp(recordingTechnique,'whole-cell') && ~strcmp(recordingTechnique,'EXCITATION') && ~strcmp(recordingTechnique,'INHIBITION')
%         if strcmp(displayName,'Led Pulse')
if strcmp(displayName,protocolID) && strcmp(egLabel,controlOrAlt) && ~strcmp(recordingTechnique,'whole-cell')
        for s = 1:length(splitFactors)
            
        splitCell{2,s}=[splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];
        end
       
        count = count + 1;
%         intensity(count) = epochs(i).meta.intensity;

        epochStorage(count,1:length(epochs(i).epoch)) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime; 
        epochStartTime(count) = epochs(i).meta.epochStartTime;
        
%         surroundSeed(count) = epochs(i).meta.surroundSeed;
   end
end

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   egLabel = epochs(i).meta.epochGroupLabel;
%    if strcmp(displayName,'Chromatic Spot') && ~strcmp(recordingTechnique,'whole-cell') && ~strcmp(recordingTechnique,'EXCITATION') && ~strcmp(recordingTechnique,'INHIBITION')
%         if strcmp(displayName,'Led Pulse')
if strcmp(displayName,protocolID) && strcmp(egLabel,controlOrAlt)
        for s = 1:length(splitFactors)
            
        splitCell{2,s}=[splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];
        end
       
        count = count + 1;
%         intensity(count) = epochs(i).meta.intensity;

        epochStorage(count,1:length(epochs(i).epoch)) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
      
        end
    end
end    

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

allSets=strings(size(splitCell{2,1},2),size(splitCell,2)-1);
clear allSets

if length(unique(splitCell{2,size(splitCell,2)})) > 1
    subtractSorter = 0;
elseif size(splitCell,2)>1
    subtractSorter = 1;
else
    subtractSorter = 0;    
end
subtractSorter=0;
for o = 1:size(splitCell,2) - subtractSorter %-1 because last parameter split is x val
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end

combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;

% if ~isnan(str2double(combos))
%     1
% combos = str2double(combos); % trying to fix the problem of sorting by first digit (strings) vs numeric value (doubles)
% combos = sort(combos);
% combos = string(combos);
% end
% while size(comboOut,1) ~= size(combos,1)


while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(comboOut,2));
    combos(1,:) = [];
end

if dataType==1
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
else

    

    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
    end
end


params = struct();

switch protocolID
case("Object Motion Texture")
    params.moveTime = 2000;
case("Object Motion Grating")
    params.moveTime = 1000;
    otherwise
    params.moveTime = 1;
end

params.saveGraph =0;
params.cellName = cellName;
params.protocolID = protocolID;
params.stimName = 'expanding';
params.dataType = dataType;

timings = [preTime stimOrig tailTime];

simpleResponse(spikeMatrix,psthMatrix,epochStorage,timings,splitCell,indexHolder,params);
%% quick plot for custom iso area summation

allAvailable = [indexHolder{1,:}];
allAvailableReshaped = reshape(allAvailable,length(indexHolder{1,1}),length(allAvailable)/length(indexHolder{1,1}));
onlyColors = allAvailableReshaped(2,:);
onlyRadii = double(allAvailableReshaped(1,:));
% if size(allAvailableReshaped,1) == 3
% onlyContrasts = double(allAvailableReshaped(3,:));
% else
%     onlyContrasts = 1;
% end

whatColors = unique(onlyColors);
whatRadii = unique(onlyRadii);
% whatContrasts = unique(onlyContrasts);

plotColors = [0 1 0;0 0 1;0 0 0];

% for a = 1:length(whatColors)
for a = 1:length(whatColors) 
    colorConditional = strcmp(whatColors(a),onlyColors);
  
    
        for b = 1:length(whatRadii)
           radiusConditional = whatRadii(b)==onlyRadii;
           finalConditional = colorConditional+radiusConditional == 2;
           if sum(finalConditional) == 0
               continue
           end
%            subplot(length(whatRadii),length(whatColors),b*a) %if for some reason more than just S- vs LM-
            figure(23)
            subplot(length(whatRadii),length(whatColors),(a+((b-1)*3))) 
           
           plot(mean(psthMatrix([indexHolder{2,finalConditional}],:))')
           title(indexHolder{1,finalConditional})

           
           spikeSumONCollection(a,b)= mean(sum(spikeMatrix([indexHolder{2,finalConditional}],timings(2)*10:(timings(2)+timings(1))*10),2));
           spikeSumONSEM(a,b) = sem(sum(spikeMatrix([indexHolder{2,finalConditional}],timings(2)*10:(timings(2)+timings(1))*10),2));
           
           spikeSumOFFCollection(a,b)= mean(sum(spikeMatrix([indexHolder{2,finalConditional}],(timings(1)+timings(2))*10:(sum(timings)*10)),2));
           spikeSumOFFSEM(a,b) = sem(sum(spikeMatrix([indexHolder{2,finalConditional}],(timings(1)+timings(2))*10:(sum(timings)*10)),2));
           holdRadii(a,b) = unique(onlyRadii(radiusConditional));
        end
        
            
       plotLength = sum(strcmp(whatColors(a),onlyColors));
%        plotLength = plotLength/2;
       
       spikeSumONCollection(spikeSumONCollection==0)=nan;
       spikeSumOFFCollection(spikeSumOFFCollection==0)=nan;

 figure(25+a)
plot(holdRadii(a,:),spikeSumONCollection(a,:),'-ok')
hold on
plot(holdRadii(a,:),spikeSumOFFCollection(a,:),'--*k')
legend('On','Off')
xlabel('spot size (micron)')
ylabel('mean spike count')
title(whatColors(a))
axis([0 inf 0 65])

        
     
     
     
     
end


%%  Expanding Spots
dataType = 1; %0 if currents
currAnalysis = 'excitation'; %exc or inh for expanding spots protocol
% currAnalysis = 'inhibition';
manualParse = false; firstTime=true; %only do manual threshold if required (turn simpleParse to false to use)

nameCurrent = "Whole cell_inh"; % either Whole cell_exc or Whole cell_inh
splitFactors = ["spotIntensity","backgroundIntensity","currentSpotSize"];

runSplitter = ["currentSpotSize"];

leakThreshold = 100;

desiredSTD =4;
% 
splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end
%%%%%%%trying to solve epoch length problem
%   counter = 0;
%   epochLength =0;
% for a = 1:length(epochs)
%    displayName = epochs(a).meta.displayName;
%    falseLength = length(epochs(a).epoch);
%    if strcmp(displayName,'Expanding Spots')
%        counter = counter +1;
%            if falseLength > epochLength
%                finale = falseLength;
%            end
%        epochLength = length(epochs(a).epoch);
%    end
% end
%%%%%%%

count = 0;
clear epochStorage
% epochStorage = zeros(counter,finale);
currType = 'no';

if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
%    egLabel = epochs(i).meta.epochGroupLabel;
%    recording = epochs(i).meta.recordingTechnique;

     if isfield(epochs(i).meta,'preTime')
      xpreTime = epochs(i).meta.preTime;
     end
  
      if isfield(epochs(i).meta,'chromaticClass')
      chroma = epochs(i).meta.chromaticClass;
     end
     
  if isfield(epochs(i).meta,'onlineAnalysis')
      oAnalysis = epochs(i).meta.onlineAnalysis;
  end
   
    if strcmp(displayName,'Expanding Spots') 
        
        
       
           oAnalysis = epochs(i).meta.onlineAnalysis;

            %solve the size problem
            
           
           leakC = epochs(i).epoch(1);
           
            %below allows manual entry of proper leak
       if manualParse   
            plot(epochs(i).epoch)
               if firstTime
                  respondedToEpoch = false; 
               else
                  respondedToEpoch = true;
               end
            while ~respondedToEpoch
                2+2
               firstTime = false;
               leakThreshold = input('enter leak: ','s');
               leakThreshold = str2num(leakThreshold);
               respondedToEpoch = true; 
            end
       end 
    end
   
   if strcmp(displayName,'Expanding Spots') && leakThreshold>leakC && leakC>-leakThreshold

       
%        if strcmp(displayName,'S MT Fspot') && strcmp(oAnalysis,'none') && xpreTime == 250 && strcmp(chroma,'blue')
       
       
        count = count + 1;
        for s = 1:length(splitFactors)
        splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
        end
        
%         width(count) = epochs(i).meta.barWidth;
        epochStorage(count,1:length(epochs(i).epoch)) = epochs(i).epoch;
%         apertureRadius(count) = epochs(i).meta.apertureRadius;
%         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
%         angleO(count) = epochs(i).meta.orientation;
   end
end

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   if strcmp(displayName,'Expanding Spots')

       leakC = epochs(i).epoch(1); 
      if manualParse
            %below allows manual entry of appropriate leak
           plot(epochs(i).epoch)
               if firstTime
                  respondedToEpoch = false; 
               else
                  respondedToEpoch = true;
               end
            while ~respondedToEpoch
               firstTime = false;
               leakThreshold = input('enter leak: ','s');
               leakThreshold = str2num(leakThreshold);
                respondedToEpoch = true; 
            end
      end
       if leakC<-leakThreshold 
           currType = 'excitation'
           leakC
       elseif leakC > leakThreshold 
           currType = 'inhibition'
       end
   end
       if strcmp(displayName,'Expanding Spots') && strcmp(currType,currAnalysis) && strcmp(nameCurrent,egLabel)

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
end
    
stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

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

if dataType == 1
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

else

    

    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
    end
end
timings = [preTime stimOrig tailTime];
params = struct();
params.saveGraph = 0;
params.stimName = 'expanding';
params.dataType = dataType;
% saveGraph = 0;
% stimName = 'Grating';

% tester = cell(2,1);
% tester{1,1}=indexHolder{1,2};
% tester{2,1}=indexHolder{2,2};
simpleSpots(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);

%% Drifting Grating
dataType =1;
% currWanted = 'excitation';
% currWanted = 'inhibition';
splitFactors = ["apertureRadius","barWidth","temporalFrequency","orientation"];
runSplitter = ["orientations"]; % not in use; logic problem

desiredSTD = 5;

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

controlOrAlt = 'Control';
% controlOrAlt = 'AltCenter';
% controlOrAlt = 'NDF 1';

count = 0;
clear epochStorage 
clear width
clear apertureRadius
clear temporalFrequency
clear angle0

if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
   end
   
   if strcmp(displayName,'Grating DSOS') && strcmp(egLabel,controlOrAlt)
%    if strcmp(displayName,'Plaid DSOS') && strcmp(egLabel,controlOrAlt)
       count=count+1;
        for s=1:length(splitCell)  
          if strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'double')
              stringedEntry = convertCharsToStrings(num2str(getfield(epochs(i).meta,splitFactors(s))));
          elseif strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'char')
              stringedEntry = convertCharsToStrings(getfield(epochs(i).meta,splitFactors(s)));
          end

            splitCell{2,s}=[splitCell{2,s} stringedEntry];
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

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis
       leakC = epochs(i).epoch(1); 
       if leakC<-50
           currType = 'excitation'
       elseif leakC > 50
           currType = 'inhibition'
       end
       
       if strcmp(egLabel,'Whole cell_exc')
           currType = 'excitation';
       elseif strcmp(egLabel,'Whole cell_inh')
           currType = 'inhibition';
       end
   end
       if strcmp(displayName,'Grating DSOS') && strcmp(oAnalysis,'analog') && strcmp(currType,currWanted)
           count = count + 1;
        for s=1:length(splitCell)  
          if strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'double')
              stringedEntry = convertCharsToStrings(num2str(getfield(epochs(i).meta,splitFactors(s))));
          elseif strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'char')
              stringedEntry = convertCharsToStrings(getfield(epochs(i).meta,splitFactors(s)));
          end

            splitCell{2,s}=[splitCell{2,s} stringedEntry];
        end
            
    %         width(count) = epochs(i).meta.barWidth;
            epochStorage(count,:) = epochs(i).epoch;
    %         apertureRadius(count) = epochs(i).meta.apertureRadius;
    %         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
%             intensity(count) = epochs(i).meta.intensity;
            preTime = epochs(i).meta.preTime;
            stimTime = epochs(i).meta.stimTime;
            tailTime = epochs(i).meta.tailTime;
    %         angleO(count) = epochs(i).meta.orientation;
       end
    end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

% for m = 1:size(splitCell,2)-1
%     if strcmp(class(splitCell{2,m}),'string')
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

%SORT THIS SHIT BABY!!!!!!
% allSets=zeros(size(splitCell{2,1},2),size(splitCell,2)-1);
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
% while ~isempty(combos)
%     counter = counter +1;
%     istherenobetterway = combos(1,:) ==allSets;
%     
%     
%     
%     indexHolder{1,counter} = combos(1,:);
%     indexHolder{2,counter} = find(sum(istherenobetterway,2)==3);
%     combos(1,:) = [];
% end

while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    
    
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(comboOut,2));
    combos(1,:) = [];
end


if dataType ==1
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

else   
    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
    end
end
timings = [preTime stimOrig tailTime];
params = struct();
params.saveGraph = 0;
params.stimName = 'Grating';
params.killCycle1 = [];  %or [] (not 0)
% timings(1)=timings(1)*2;
params.time = timings;
% params.killCycle1 = 'false';
% saveGraph = 0;
% stimName = 'Grating';

orientedStim(epochStorage,spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

%% Oriented Bars
dataType = 1; %0 if currents
currAnalysis = 'excitation';
% currAnalysis = 'inhibition';
splitFactors = ["intensity","backgroundIntensity","barSize","orientation"];
% splitFactors = ["contrast","backgroundIntensity","chromaticClass","barSize","orientation"];

%NOTE: last split is always X axis.  I think this is helpful because can be
%specified by length function and don't need to cary another thing to
%processing function
subtractBGrate=1;

desiredSTD =4;
% 
protocolType = 'Oriented Bars';
% protocolType = 'Oriented Bars Fixed';
% protocolType = 'Chromatic Bars';

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

controlOrAlt = "Control";
% controlOrAlt = "NDF 4";  
% controlOrAlt =  "AltCenter";

% egLabelforWholeCell = "Wholecell_exc";

count = 0;
clear epochStorage;
count2 =0;
if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
     egLabel = epochs(i).meta.epochGroupLabel;
    if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
    end
    
    if isfield(epochs(i).meta,'centerOffsetTrue')
       offset = epochs(i).meta.centerOffsetTrue;
    end
    
   
              if strcmp(displayName,protocolType) && strcmp(egLabel,controlOrAlt)
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
            preTime = epochs(i).meta.preTime;
            stimTime = epochs(i).meta.stimTime;
            tailTime = epochs(i).meta.tailTime;
              end

              
end
       

        
      
 


else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recordingTechnique = epochs(i).meta.recordingTechnique;
      egLabel = epochs(i).meta.epochGroupLabel
%    if isfield(epochs(i).meta,'onlineAnalysis')
%        oAnalysis = epochs(i).meta.onlineAnalysis
%        leakC = epochs(i).epoch(1)
%        if leakC<-5
%            currType = 'excitation'
%        elseif leakC > 0
%            currType = 'inhibition'
%        end
%    end
%    & strcmp(currType,currAnalysis)
       if strcmp(displayName,protocolType) && strcmp(egLabel,egLabelforWholeCell) 
      
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
            preTime = epochs(i).meta.preTime;
            stimTime = epochs(i).meta.stimTime;
            tailTime = epochs(i).meta.tailTime;
       end
    end
end

stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

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

if dataType ==1
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

else   
    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:5000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:5000));
    end
end

params = struct();
params.saveGraph =1;
params.saveIter = 1;
params.stimName = 'bars';
params.dataType = dataType;
params.bgRate = subtractBGrate;
params.fileName = cellName;
timings = [preTime stimOrig tailTime];
orientedStim(epochStorage,spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% MOVING BAR

dataType = 1; %0 if currents
splitFactors = ["intensity","backgroundIntensity","speed","barSize","orientation"]; colorBar = false;
% splitFactors = ["barColor","barSize","speed","orientation"]; colorBar = true;
%NOTE: last split is always X axis.  I think this is helpful because can be
%specified by length function and don't need to cary another thing to
%processing function




leakC = 100;
splitCell = cell(2,length(splitFactors));

desiredSTD = 4;

if colorBar
    desiredProtocol = 'Moving Chromatic Bar';
else
    desiredProtocol = 'Moving Bar';
end


% controlOrAlt = "Wholecell_exc";
controlOrAlt = "Control";

% controlOrAlt = "AltCenter";


for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end


count = 0;
clear epochStorage;
if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
    egLabel = epochs(i).meta.epochGroupLabel;
  if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
   end
   
%    if strcmp(displayName,'Moving Bar') && strcmp(egLabel,'Control') && strcmp(oAnalysis,'extracellular')



if strcmp(displayName,desiredProtocol) && strcmp(egLabel,controlOrAlt)
                for s=1:length(splitCell)  
                  if strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'double')
                      stringedEntry = convertCharsToStrings(num2str(getfield(epochs(i).meta,splitFactors(s))));
                  elseif strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'char')
                      stringedEntry = convertCharsToStrings(getfield(epochs(i).meta,splitFactors(s)));
                  end
               
                    splitCell{2,s}=[splitCell{2,s} stringedEntry];
                end
       
        count = count + 1;
%         intensity(count) = epochs(i).meta.intensity;
        epochStorage(count,:) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
      
   end
end

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
       leakC = epochs(i).epoch(1); 
   end

       if strcmp(displayName,desiredProtocol) && strcmp(egLabel,controlOrAlt)
%            if strcmp(displayName,'Moving Bar') && ~strcmp(oAnalysis,'extracellular')
%            leakC
            barSizeFind = contains(splitFactors,'barSize');
                for s=1:length(splitCell)  
                  if strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'double')
                      stringedEntry = convertCharsToStrings(num2str(getfield(epochs(i).meta,splitFactors(s))));
                  elseif strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'char')
                      stringedEntry = convertCharsToStrings(getfield(epochs(i).meta,splitFactors(s)));
                  end
               
                    splitCell{2,s}=[splitCell{2,s} stringedEntry];
                end
            count = count + 1;
    %         width(count) = epochs(i).meta.barWidth;
            epochStorage(count,:) = epochs(i).epoch;
    %         apertureRadius(count) = epochs(i).meta.apertureRadius;
    %         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
%             intensity(count) = epochs(i).meta.intensity;
            preTime = epochs(i).meta.preTime;
            stimTime = epochs(i).meta.stimTime;
            tailTime = epochs(i).meta.tailTime;
    %         angleO(count) = epochs(i).meta.orientation;
       end
    end
end




stimOrig = stimTime;
sampleRate = 10000;
% preTime = epochs(1).meta.preTime;
% stimTime = epochs(1).meta.stimTime;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));

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

if dataType == 1
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

else   
    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:5000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:5000));
    end
end

params = struct();
params.saveGraph =0;
%grating,bars
params.stimName = 'moving bar';
params.saveIter = 1;
params.Region = 7500:11000;
% params.Region = 13000:17000;
% params.Region = 1:35000;
timings = [preTime stimOrig tailTime];
dataToAlign = orientedStim(epochStorage,spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% plottin' code

            figure(10); hold on
           
            for k = 1:12
            plot((dataToAlign(k,:)/max(abs(dataToAlign(k,:))))+(1*(k-1)),'b')
            plot((excData(k,:)/max(abs(excData(k,:))))+(1*(k-1)),'r') 
            end
title('Moving Bar: Normalized Spiking PSTHs & Normalized Excitation at Each Angle')
                        figure(11); hold on
           
            for k = 1:12
            plot(excData(k,:)+(1000*(k-1)),'b','LineWidth',1)
            plot(inhibData(k,:)+(1000*(k-1)),'r','LineWidth',2) 
            end
legend('excitation','inhibition')
title('Moving Bar: Inhibition & Excitation at Each Angle')
 xlabel('time pts (2.5 second bar presentation)')
set(gca,'ytick',[])
%% processing functions
    
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

function CRF(spikeMatrix,psthMatrix,radius,timings,saveStuff)


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

% if saveStuff.flag == 1 && strcmp(saveStuff.cellType,'OFF Parasol')
%     disp('accessed')
%     MNOFFParasol.Exp(saveStuff.expNum).Cell(saveStuff.OFFParasolC).Sensitivity = [xaxis;yaxis];
% elseif saveStuff.flag == 1 && strcmp(saveStuff.cellType,'ON Parasol')
%     MNONParasol.Exp(saveStuff.expNum).Cell(saveStuff.ONParasolC).Sensitivity = [xaxis;yaxis];
% end
plot(xaxis,yaxis)
axis([0 max(xaxis) 0 max(yaxis)+2])

end



