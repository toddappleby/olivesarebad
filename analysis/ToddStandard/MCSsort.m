%% MCS

dataType = 1; %0 if currents
splitFactors = ["centerBarOrientation","surroundBarOrientation","surroundBarFrameDwell","centerBarWidth","contrast"];
%NOTE: last split is always X axis.  I think this is helpful because can be
%specified by length function and don't need to cary another thing to
%processing function
leakC = 0;
splitCell = cell(2,length(splitFactors));

desiredSTD = 6;


for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

%this should solve epochs of different lengtsh
  counter = 0;
  epochLength =0;
for a = 1:length(epochs)
   displayName = epochs(a).meta.displayName;
   falseLength = length(epochs(a).epoch);
   if strcmp(displayName,'Motion Center Surround')
       counter = counter +1;
           if falseLength > epochLength
               finale = falseLength;
           end
       epochLength = length(epochs(a).epoch);
   end
end

epochStorage = zeros(counter,finale);

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
   
   if strcmp(displayName,'Motion Center Surround') && strcmp(oAnalysis,'extracellular') && ~strcmp(egLabel,'Bad Focus')
                for s=1:length(splitCell)  
                  if strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'double')
                      stringedEntry = convertCharsToStrings(num2str(getfield(epochs(i).meta,splitFactors(s))));
                  elseif strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'char')
                      stringedEntry = convertCharsToStrings(getfield(epochs(i).meta,splitFactors(s)));
                  end
               
                    splitCell{2,s}=[splitCell{2,s} stringedEntry];
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

       if strcmp(displayName,'Motion Center Surround') && strcmp(oAnalysis,'analog') && leakC > -500
           
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
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
    end
end

params = struct();
params.saveGraph =0;
params.stimName = 'MCS';
timings = [preTime stimOrig tailTime];
MCS(epochStorage,spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
