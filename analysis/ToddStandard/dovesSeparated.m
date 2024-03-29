% load data -- #2
%  cd('E:\Data Analysis_2020\2020_0611\')
% cd('C:\Users\reals\Documents\PhD 2021\ClarinetExports\2019_1022')
%  fileName = '20191022Bc2.mat';
% saveFile = "20191022Bc2_doves";
 fileName = '20221219Ac1.mat';
saveFile = "20221219Ac1_doves";
 load(fileName)

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

%% Doves
dataType = 1; %0 if currents
nameCurrent = "Whole cell_exc"; % either Whole cell_exc or Whole cell_inh
splitFactors = ["maskDiameter","apertureDiameter","stimulusIndex"];
desiredSTD = 6;
splitCell = cell(2,length(splitFactors));

for g = 1:length(splitCell)
    splitCell{1,g} = splitFactors(g);
end

count = 0;
clear epochStorage 

controlOrAlt = 'Control';


if dataType == 1
    subjectName = {};
    imageName = {};
for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
   
    if strcmp(displayName,'Doves Movie') && strcmp(egLabel,controlOrAlt)
%        oAnalysis = epochs(i).meta.onlineAnalysis;
       leakC = epochs(i).epoch(1); 
       if 50>leakC && leakC>-50
           currType = 'yes'
           leakC
       else
           currType = 'no'
           leakC
       end
    end
   currType='yes';
   if strcmp(displayName,'Doves Movie') && strcmp(currType,'yes')
        
        count = count + 1;
        for s = 1:length(splitFactors)
        splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
        end
        
        %data
        epochStorage(count,:) = epochs(i).epoch;
        
        %timings
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
        
        %extras
        subjectName{count} = epochs(i).meta.subjectName;
        imageName{count} = epochs(i).meta.stimulusIndex;
       
   end
end

else
    for i = 1:length(epochs)
  
   displayName = epochs(i).meta.displayName;
   egLabel = epochs(i).meta.epochGroupLabel;
   recording = epochs(i).meta.recordingTechnique;
       if strcmp(displayName,'Doves Movie') && strcmp(egLabel,nameCurrent)

            count = count + 1;
            for s = 1:length(splitFactors)
            splitCell{2,s}=[splitCell{2,s} getfield(epochs(i).meta,splitFactors(s))];
            end

        %data
        epochStorage(count,:) = epochs(i).epoch;
        
        %timings
        preTime = epochs(i).meta.preTime;
        stimTime = epochs(i).meta.stimTime;
        tailTime = epochs(i).meta.tailTime;
        
        %extras
        subjectName{count} = epochs(i).meta.subjectName;
        imageName{count} = epochs(i).meta.imageName;
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
            psthMatrix(k,:) = psth(spikeMatrix(k,:),10,sampleRate,1);
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
params.saveFile = saveFile;
% saveGraph = 0;
% stimName = 'Grating';
DOVES(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params);
%%
function DOVES(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

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
            plot(psthData(b,:)+(240*(b-1)))
            title(indexHolder{1,a})
            figure(11); hold on
            subplot(1,size(indexHolder,2),a)
            spikeData(b,:) = sum(spikeMatrix(sortedIndex(finalInd),:),1);
            plot(spikeData(b,:)+(2*(b-1)))
            title(indexHolder{1,a})
            else
            psthData(b,:) = mean(spikeMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(b,:)+(100*(b-1)))
            title(indexHolder{1,a})
            end
            
            
            
         end
        
         
         

      
    end
    
    uniqueIndex = unique(splitCell{2,size(splitCell,2)});
    stimIndex = {};
%     load('E:\Data Analysis_2020\code\Manookin Repository\manookin-package\resources\dovesFEMstims20160826.mat')
     load('C:\Users\reals\Documents\manookin-package\resources\doves\dovesFEMstims20160826.mat')
labels = ["center","fullfield","surround"];
cellType = "RBS_";

% for e = 1:length(indexHolder)
for e = 1:3
 
    for g = 1:length(uniqueIndex)
        
        picIndex = find(splitCell{2,size(splitCell,2)}==uniqueIndex(g));

        stimPlacer = ismember(indexHolder{2,1},picIndex); % not sure how to automate this part, hopefully remains the same each time (to get full field image only);
        stimIndex{1,g} = uniqueIndex(g);
        stimIndex{2,g} = FEMdata(uniqueIndex(g)).ImageName;
        stimIndex{3,g} = mean(psthMatrix(indexHolder{2,1}(stimPlacer),:),1);
%         forIgor(g,:) = mean(psthMatrix(indexHolder{2,e}(stimPlacer),:),1); %e for center/full/surround Note: it just happens to fall in this order -- not automated
%         size(forIgor)
    end
    saveName = strcat('C:\Users\reals\Documents\PhD 2021\ClarinetExports\',params.saveFile,cellType,labels(e),".mat");
%       saveName = strcat('E:\Data Analysis_2020\weeklymeeting_529\',params.saveFile,cellType,labels(e),".mat");
%      save(saveName,'forIgor')
    
end
save(saveName,'stimIndex')



end
