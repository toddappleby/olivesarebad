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
cd('E:\Data Analysis_2020\2019_0620\')
expDate = dir;
load('20190620Bc3.mat')
frameTimings = epochs;
%% what protocols contained in data?

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

%% chirp regurgitation
dataType = 1;
clear epochStorage

count =0;

 
if dataType == 1
for i = 1:length(epochs)
   displayName = epochs(i).meta.displayName;

   recordingTechnique = epochs(i).meta.recordingTechnique;

   egLabel = epochs(i).meta.epochGroupLabel;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
   end


   if strcmp(displayName,'Chirp Stimulus') && ~strcmp(recordingTechnique,'whole-cell') && strcmp(egLabel,'Control')
count = count + 1;
epochStorage(count,:) = epochs(i).epoch;

% contrast(count) = epochs(i).meta.contrast;

preTime = epochs(i).meta.preTime;

stimTime = epochs(i).meta.stimTime;

tailTime = epochs(i).meta.tailTime;

   

 

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

preTime = epochs(i).meta.preTime;

stimTime = epochs(i).meta.stimTime;

tailTime = epochs(i).meta.tailTime;

   

 

        end

    end  
end

for z = 1:size(epochStorage,1)

plot(epochStorage(z,:))

hold on

end


    

%% MTF spots and annuli

splitFactors = ["stimulusClass","temporalFrequency","temporalClass","radius"];
runSplitter = ["radius"]; 

splitCell = cell(2,length(splitFactors));

desiredSTD = 6;

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
   
  if isfield(epochs(i).meta,'onlineAnalysis')
      oAnalysis = epochs(i).meta.onlineAnalysis;
  end
   
   if strcmp(displayName,'S MT Fspot') && strcmp(oAnalysis,'extracellular') 
    temporalClass = epochs(i).meta.temporalClass;
    if strcmp(temporalClass,'sinewave')
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

%DON'T GIVE A FUCK ABOUT A STRING!!!!

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
desiredSTD = 8;
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
dataType = 1; %0 if currents
currAnalysis = 'excitation'; %exc or inh for expanding spots protocol
% currAnalysis = 'inhibition';
nameCurrent = "Whole cell_inh"; % either Whole cell_exc or Whole cell_inh
splitFactors = ["spotIntensity","backgroundIntensity","currentSpotSize"];
runSplitter = ["currentSpotSize"]; %might just set this manually because not gonna change? number could change if the protocol run was cut short
%this is why Greg did this the way he did .... maybe.  gonna have to figure
%this out by TIME??? which seems rough

desiredSTD = 5;

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage 
currType = 'no';
if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
    if strcmp(displayName,'Expanding Spots')
       oAnalysis = epochs(i).meta.onlineAnalysis;
       leakC = epochs(i).epoch(1); 
       if 50>leakC && leakC>-50
           currType = 'yes'
           leakC
       else
           currType = 'no'
       end
   end
   
   if strcmp(displayName,'Expanding Spots') && strcmp(currType,'yes')
        
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

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
   if strcmp(displayName,'Expanding Spots')
%        oAnalysis = epochs(i).meta.onlineAnalysis
       leakC = epochs(i).epoch(1); 
       if leakC<-50
           currType = 'excitation'
       elseif leakC > 50
           currType = 'inhibition'
       end
   end
       if strcmp(displayName,'Expanding Spots') && strcmp(currType,currAnalysis)

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
simpleSpots(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);

%% Drifting Grating
dataType =1;
currWanted = 'excitation';
% currWanted = 'inhibition';
splitFactors = ["apertureRadius","barWidth","temporalFrequency","orientation"];
runSplitter = ["orientations"]; %might just set this manually because not gonna change? number could change if the protocol run was cut short
%this is why Greg did this the way he did .... maybe.  gonna have to figure
%this out by TIME??? which seems rough

desiredSTD = 5;


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

if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
   end
   
   if strcmp(displayName,'Grating DSOS') && strcmp(oAnalysis,'extracellular')    
        
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

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis
       leakC = epochs(i).epoch(1); 
       if leakC<-50
           currType = 'excitation'
       elseif leakC > 50
           currType = 'inhibition'
       end
   end
       if strcmp(displayName,'Grating DSOS') && strcmp(oAnalysis,'analog') && strcmp(currType,currWanted)
           
            barSizeFind = contains(splitFactors,'barSize');
          for s = 1:length(splitFactors)
            if barSizeFind(s)
                count2=count2+1
               bSizeDisco = getfield(epochs(i).meta,splitFactors(s)); %disonnected bar size (as length 2 array) 
               strTransferTicket = strcat(num2str(bSizeDisco(1)),num2str(bSizeDisco(2)))
               barSizeCombined = str2num(strTransferTicket);
               splitCell{2,s}=[splitCell{2,s} barSizeCombined];
            else
            splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
           
            end
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
% saveGraph = 0;
% stimName = 'Grating';
orientedStim(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% Oriented Bars
dataType = 1; %0 if currents
currAnalysis = 'excitation';
% currAnalysis = 'inhibition';
splitFactors = ["intensity","backgroundIntensity","barSize","orientation"];
%NOTE: last split is always X axis.  I think this is helpful because can be
%specified by length function and don't need to cary another thing to
%processing function

desiredSTD = 4;

splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end


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
   
   if strcmp(displayName,'Oriented Bars') && strcmp(oAnalysis,'extracellular')    
        
        for s = 1:length(splitFactors)
            barSizeFind = contains(splitFactors,'barSize');
            if barSizeFind(s)
                count2=count2+1;
               bSizeDisco = getfield(epochs(i).meta,splitFactors(s)); %disonnected bar size (as length 2 array) 
               strTransferTicket = strcat(num2str(bSizeDisco(1)),num2str(bSizeDisco(2)))
               barSizeCombined = str2num(strTransferTicket);
               splitCell{2,s}=[splitCell{2,s} barSizeCombined];
            else
               splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
            end
        
        end
       
        count = count + 1;
        intensity(count) = epochs(i).meta.intensity;
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
   recording = epochs(i).meta.recordingTechnique;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis
       leakC = epochs(i).epoch(1); 
       if leakC<-50
           currType = 'excitation'
       elseif leakC > 50
           currType = 'inhibition'
       end
   end
       if strcmp(displayName,'Oriented Bars') && strcmp(oAnalysis,'analog') && strcmp(currType,currAnalysis)
            barSizeFind = contains(splitFactors,'barSize');
          for s = 1:length(splitFactors)
            if barSizeFind(s)
                count2=count2+1
               bSizeDisco = getfield(epochs(i).meta,splitFactors(s)); %disonnected bar size (as length 2 array) 
               strTransferTicket = strcat(num2str(bSizeDisco(1)),num2str(bSizeDisco(2)))
               barSizeCombined = str2num(strTransferTicket);
               splitCell{2,s}=[splitCell{2,s} barSizeCombined];
            else
            splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
           
            end
          end
            count = count + 1;
    %         width(count) = epochs(i).meta.barWidth;
            epochStorage(count,:) = epochs(i).epoch;
    %         apertureRadius(count) = epochs(i).meta.apertureRadius;
    %         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
            intensity(count) = epochs(i).meta.intensity;
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

params = struct();
params.saveGraph =0;
params.stimName = 'bars';
timings = [preTime stimOrig tailTime];
orientedStim(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%% MOVING BAR

dataType = 1; %0 if currents
splitFactors = ["intensity","backgroundIntensity","barSize","orientation"];
%NOTE: last split is always X axis.  I think this is helpful because can be
%specified by length function and don't need to cary another thing to
%processing function
leakC = 0;
splitCell = cell(2,length(splitFactors));

desiredSTD = 5;


for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage;
if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
  if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis;
   end
   
   if strcmp(displayName,'Moving Bar') && strcmp(oAnalysis,'extracellular')      
        for s = 1:length(splitFactors)
            barSizeFind = contains(splitFactors,'barSize');
            if barSizeFind(s)
               bSizeDisco = getfield(epochs(i).meta,splitFactors(s)); %disonnected bar size (as length 2 array) 
               strTransferTicket = strcat(num2str(bSizeDisco(1)),num2str(bSizeDisco(2)));
               barSizeCombined = str2num(strTransferTicket);
               splitCell{2,s}=[splitCell{2,s} barSizeCombined];
            else
               splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
            end
        end
       
        count = count + 1;
        intensity(count) = epochs(i).meta.intensity;
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
   recording = epochs(i).meta.recordingTechnique;
   if isfield(epochs(i).meta,'onlineAnalysis')
       oAnalysis = epochs(i).meta.onlineAnalysis
       leakC = epochs(i).epoch(1); 
   end
       if strcmp(displayName,'Moving Bar') && strcmp(oAnalysis,'analog') && leakC < -50
           
            barSizeFind = contains(splitFactors,'barSize');
          for s = 1:length(splitFactors)
            if barSizeFind(s)
                count2=count2+1
               bSizeDisco = getfield(epochs(i).meta,splitFactors(s)); %disonnected bar size (as length 2 array) 
               strTransferTicket = strcat(num2str(bSizeDisco(1)),num2str(bSizeDisco(2)))
               barSizeCombined = str2num(strTransferTicket);
               splitCell{2,s}=[splitCell{2,s} barSizeCombined];
            else
            splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
           
            end
          end
            count = count + 1;
    %         width(count) = epochs(i).meta.barWidth;
            epochStorage(count,:) = epochs(i).epoch;
    %         apertureRadius(count) = epochs(i).meta.apertureRadius;
    %         temporalFrequency(count) = epochs(i).meta.temporalFrequency;
            intensity(count) = epochs(i).meta.intensity;
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
            
            if params.dataType ==1
            psthData(b,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(b,:)+(100*(b-1)))
            title(indexHolder{1,a})
            else
            psthData(b,:) = mean(spikeMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(b,:)+(100*(b-1)))
            title(indexHolder{1,a})
            end
            
            
            
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