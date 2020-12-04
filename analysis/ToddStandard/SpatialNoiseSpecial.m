%% Spatial  Noise

splitFactors = ["frameDwell","noiseClass"];


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
   
              if strcmp(displayName,'Spatial Noise') && strcmp(oAnalysis,'extracellular')  
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
            numXChecks = epochs(i).meta.numXChecks;
            numYChecks = epochs(i).meta.numYChecks;
            seeds(count)= epochs(i).meta.seed;
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
        psthMatrix(k,:) = psth(spikeMatrix(k,:),6+2/3,sampleRate,10);
          end
end

%60 f/sec * 21sec 
numFrames = 1260;

params = struct();
params.saveGraph =0;
params.stimName = 'bars';
timings = [preTime stimOrig tailTime];
[allM,frameVals] = doSpatialMap(psthMatrix,numXChecks,numYChecks,numFrames,seeds,timings);
%% plot test 1

for c = 1:size(allM,3)
    imagesc(allM(:,:,c))
    pause
end

%% plot test 2
for p = 1:6
subplot(2,3,p)
imagesc(allM(:,:,p+4))
title(p+1)
end
