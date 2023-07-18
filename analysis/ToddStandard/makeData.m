function [splitCell,indexHolder,spikeMatrix,frameTimings,meta,seed] = makeData(epochs,frameTs,protocolID,splitFactors,desiredSTD,recording,newDataFlag)
%%%%setup
if strcmp(recording,'extracellular')
    dataType = 1; 
else
    dataType = 0;
end

splitCell = cell(2,length(splitFactors));

for g = 1:size(splitCell,2)
    splitCell{1,g} = splitFactors(g);
end

  counter = 0;
  epochLength =0;
  
for a = 1:length(epochs)
   displayName = epochs(a).meta.displayName;
   falseLength = length(epochs(a).epoch);
   if strcmp(displayName,protocolID)
       counter = counter +1;
       
           if falseLength > epochLength
               finale = falseLength;
           end
           
       epochLength = length(epochs(a).epoch);
   end
end

epochStorage = zeros(counter,finale);
count = 0;

%%%%% start
if dataType == 1
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   egLabel = epochs(i).meta.epochGroupLabel;
     if strcmp(egLabel,'motion and noise')
         egCompare = 'motion and noise';
     elseif strcmp(egLabel, 'AltCenter')
         egCompare = 'AltCener';
     else
         egCompare = 'Control';
   end
% egCompare = 'AltCenter';
if strcmp(displayName,protocolID) && strcmp(egLabel,egCompare) && ~strcmp(recordingTechnique,'whole-cell')
        for s = 1:length(splitFactors)
   
            if strcmp(splitFactors(s),'frameDwell') %correction for M&N where framedwell wasn't always a parameter
                 if isfield(epochs(i).meta,'frameDwell') 
                     splitCell{2,s}=[splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];
                 else
                     splitCell{2,s}=[splitCell{2,s} "1"];
                 end
            else
                if strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'double')
                      stringedEntry = convertCharsToStrings(num2str(getfield(epochs(i).meta,splitFactors(s))));
                  elseif strcmp(class(getfield(epochs(i).meta,splitFactors(s))),'char')
                      stringedEntry = convertCharsToStrings(getfield(epochs(i).meta,splitFactors(s)));
                end
                  splitCell{2,s}=[splitCell{2,s} stringedEntry];
%                 splitCell{2,s}=[splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];
            end
             
        
        end
        
        
       
        count = count + 1;
        if strcmp(protocolID,'Motion And Noise')
        seed(count) = epochs(i).meta.seed;
        end
        epochStorage(count,1:length(epochs(i).epoch)) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime; 
        holdi = i;
      
   end
end

meta = epochs(holdi).meta;

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   recordingTechnique = epochs(i).meta.recordingTechnique;
   egLabel = epochs(i).meta.epochGroupLabel;

if strcmp(displayName,protocolID) && strcmp(egLabel,'Whole cell_exc')
        
    
    for s = 1:length(splitFactors)
            
            splitCell{2,s}=[splitCell{2,s} string(getfield(epochs(i).meta,splitFactors(s)))];
        
        end
       
        count = count + 1;

        epochStorage(count,1:length(epochs(i).epoch)) = epochs(i).epoch;
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
      
        end
    end
end    

stimOrig = stimTime;
sampleRate = 10000;
stimTime = [preTime/(1/10000) (preTime+stimTime)/(1/10000)];
stimTime = stimTime/1000;

spikeMatrix = zeros(size(epochStorage,1),size(epochStorage,2));
psthMatrix = zeros(size(epochStorage,1),size(epochStorage,2));


%%%%start of big sort
allSets=strings(size(splitCell{2,1},2),size(splitCell,2)-1);
clear allSets


if size(splitCell,2)>1
    subtractSorter = 1;
else
    subtractSorter = 0;
end

subtractSorter = 0;

for o = 1:size(splitCell,2) - subtractSorter %-1 because last parameter split is x val
    for p = 1:size(splitCell{2,1},2)
        allSets(p,o) = splitCell{2,o}(1,p);       
    end
end

combos = unique(allSets,'rows');
comboOut = combos;
clear indexHolder
counter = 0;


while ~isempty(combos)
    counter = counter +1;
    istherenobetterway = combos(1,:) ==allSets;
    indexHolder{1,counter} = combos(1,:);
    indexHolder{2,counter} = find(sum(istherenobetterway,2)==size(comboOut,2));
    combos(1,:) = [];
end
frameTimings= [];


for f = 1:length(frameTs)

     displayName = frameTs(f).meta.displayName;
    
     
   
   if strcmp(displayName,protocolID)
       
    count = count+1;
    singleMonitorRun = frameTs(f).epoch;
 frameTimings(count,:) = singleMonitorRun;

%  epochStartTime(count,:) = frameTimings(f).meta.epochStartTime;
%  epochNumFrames(count) = frameTimings(f).meta.epochNum;
   end
end

if dataType==1 && newDataFlag==1
    %%% spike detection
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

elseif dataType == 0 %here so data flag = 0 doesn't jump to current collection    

    for k = 1:size(epochStorage,1)
        spikeMatrix(k, :) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
        psthMatrix(k,:) = epochStorage(k, :) - mean(epochStorage(k, 1:1000));
    end
    
end

if ~strcmp(protocolID,'Motion And Noise')
seed =[];
end

end



