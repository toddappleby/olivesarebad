%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'S MT Fspot';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[]; %removes macOS hidden stuff
%%
 timings = [500,2500,500]; %same for these
 sampleRate = 10000;

stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;
 nonlinearityBins = 100;
kmeansGroups=[];
cellnameList = [];
errorGroups = [];
recordingType = 'extracellular';
binRate = 1e3;
yCrossCells=[];
holdName=[];
holdFile=[];
f1collection=[];
cellF1F2=[];
cellF1=[];
cellF2=[];
phaser=[];
count=0;
areaCurve=[];
areaRadii=[];
cellCollect = cell(size(folderSet,3),2);


for c = 1:size(folderSet,3)
    
    cd(strcat(typesFolder,folderSet(:,:,c)))
    
    cellType = strtrim(folderSet(:,:,c))
    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    cellTypeData = [];
    
holdY=[];
collectedMeans=[];
collectedRadii=[];

        for d = 1:size(fileIndex,3)
      count=count+1;

         load(strtrim(fileIndex(:,:,d)))
         %Where params are most appropriate?  prob the ones I ran the most:
         
             timings = [500,500,500]; %the metadata isn't saved to smtf mats for some reason and I don't feel like fixing it right now (midnight 1/20/2022)
             if size(spikingData,2)<10001
                 timings=[250,500,250];
             end
         stimStart = (timings(1)*1e-3)*sampleRate+1;
         stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
         stimOff = (timings(1)+timings(2)+timings(3))*10;
         
         sumsort=[];
         clear labelHolder
 for s = 1:size(indexHolder,2)
     labelHolder(s,:) = string(indexHolder{1,s});
 end
 %%%%%%%% Pulse Processing
   
%        
%        
%              
 
         spotIntensityOnly = double(labelHolder(:,1));
         stimType = labelHolder(:,2);
         temporalType = labelHolder(:,3);
         posContrast = spotIntensityOnly>0;
         negContrast = spotIntensityOnly<0;
        
         spotOnly = strcmp(stimType,"spot");
         pulseOnly = strcmp(temporalType,"pulse");

         
         indPos= vertcat(indexHolder{2,spotOnly&pulseOnly&posContrast});
         indNeg = vertcat(indexHolder{2,spotOnly&pulseOnly&negContrast});
         
         if isempty(indPos) && isempty(indNeg)
             continue
         end
         
         indPoswithBGRadii = double(splitCell{2,4}(indPos));
         indNegwithBGRadii = double(splitCell{2,4}(indNeg));
        
         if length(indPoswithBGRadii) > length(indNegwithBGRadii)
             radii=unique(indPoswithBGRadii);
         else
             radii=unique(indNegwithBGRadii);
         end
         
%          radii=unique(double(splitCell{2,4}));
         radiiMean = zeros(size(radii));
         
         for xx = 1:length(radii)
            
            if length(indPos)>length(indNeg)
            arrayIndex= indPos(indPoswithBGRadii==(radii(xx)),1);
            radiiMean(xx) = mean(sum(spikingData(arrayIndex,stimStart:stimOff),2));
            else
            arrayIndex= indNeg(indNegwithBGRadii==(radii(xx)),1);
            radiiMean(xx) = mean(sum(spikingData(arrayIndex,stimStart:stimOff),2));
            end

         end
         
         if sum(radiiMean)==0 
             count = count -1;
             continue   
         end
         
         if sum(isnan(radiiMean))>0
             count=count-1;
             continue
         end
             
         
         afterMax = radiiMean(find(radiiMean==max(radiiMean)):end);
         afterMaxRadii = radii(find(radiiMean==max(radiiMean)):end);
         percentReduction = 1 - min(afterMax)/max(radiiMean);
          suppressTest = mean(radiiMean((length(radiiMean)-2):end));
         suppressionIndex = (max(radiiMean)-suppressTest)/max(radiiMean);
         maxRadius = min(radii(radiiMean==max(radiiMean)));
         
         halfMaxRadius = radii(round(find(radiiMean==max(radiiMean))/2));
         riseRateLate = (radiiMean(find(radiiMean==max(radiiMean))) - radiiMean(round(find(radiiMean==max(radiiMean))/2)))/ (maxRadius-halfMaxRadius);
         riseRateEarly = (radiiMean(round(find(radiiMean==max(radiiMean))/2)) - radiiMean(1)) / (halfMaxRadius-radii(1));
         suppressRate = (radiiMean(find(radiiMean==max(radiiMean))) - min(afterMax))/ (maxRadius - min(afterMaxRadii(afterMax == min(afterMax))));
         
         
         sMTFSpotsOut{1,count} = string(cellType);
         sMTFSpotsOut{2,count} = string(fileIndex(:,1:11,d));
         sMTFSpotsOut{3,count} = [maxRadius(1), halfMaxRadius(1), percentReduction(1), riseRateEarly(1), riseRateLate(1), suppressRate(1) suppressionIndex(1)];
         
          collectedMeans(d,1:length(radiiMean)) = radiiMean;
          collectedRadii(d,1:length(radii))=radii;
          
          
          
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
          %SinewaveProcessing
%     sumsort = sum(contains(labelHolder,["1" "annulus" "sinewave"]),2);
%        if sumsort~=2
%            continue
%        end
%        indicesWanted = indexHolder{2,sumsort==2};
%         
%         radii = splitCell{2,4};
%         radii = str2double(radii); 
%         indexedRadii=radii(indicesWanted);
%         
%         avgF1=[];
%         avgF2=[];
%         
%       if size(spikingData,2)<30000
%           continue
%       end
% 
%    
%         for w = 1:length(unique(indexedRadii))
%             listRadii = unique(indexedRadii);
%             data = mean(spikingData(indicesWanted(indexedRadii==listRadii(w)),:),1);
%                    binnedData = BinSpikeRate(data(stimStart:stimEnd), 100, sampleRate);
%                 [F, phase] = frequencyModulation(binnedData, ...
%                 100, 2, 'avg', 1:2, []);
%             
%                 avgF1(w)= F(1);
%                 avgF2(w)= F(2);   
%                 phaseBit(w)=phase(1);
%         end
% 
%         if length(avgF1)~=18
%             continue
%         end
% 
%       
%         if isnan(avgF1)
%             continue
%         end
% 
%         
%         avgF1 = avgF1/max(avgF1);
%         
%         cellF1 = [cellF1; avgF1];
%         cellF2 = [cellF2; avgF2];
%         phaser = [phaser; phaseBit(1:size(phaser,2))];
%         cellF1F2 = [cellF1F2; [avgF1 avgF2]];
%         cellradii=listRadii;
%         holdName= [string(holdName); string(cellType)];
%         holdFile = [string(holdFile); fileIndex(:,1:11,d)];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%  

  
        end
        onlyData=[];
        onlyDataRadii=[];
        uniqueRadii=[];
        radiiDataMeans = [];
        onlyData = collectedMeans(sum(collectedMeans,2)~=0,:);
        onlyData = onlyData./max(onlyData')';
        onlyDataRadii = collectedRadii(sum(collectedMeans,2)~=0,:);
        uniqueRadii = unique(onlyDataRadii);
        
        for orderRadii = 1:length(uniqueRadii)
               radiiDataMeans(orderRadii) =  mean(onlyData(onlyDataRadii == uniqueRadii(orderRadii)));
        end
        
        cellCollect{c,1} = radiiDataMeans;
        cellCollect{c,2} = uniqueRadii';
        cellCollect{c,3} = string(strtrim(cellType));
        
           meanType = mean(collectedMeans,1);
        areaCurve(c,1:length(meanType))=meanType;
        
        %must clear zeros from above because matlab puts zeros to
        %compensate for changing size of array #thanksmatlab (actually)
%         areaCurve(c,length(meanType):size(areaCurve,2))=areaCurve(c,length(meanType)); 
        
        areaRadii(c,1:length(radii))=radii;
        areaRadii(c,areaRadii(c,:)==0)=max(areaRadii(c,:));
       
         
         
        holdName = [string(holdName); string(cellType)];
        
        
        
end
% 
%  cd(typesFolder)
% 
% cellF1F2 = normr(cellF1F2); 
% cellF1 = normr(cellF1); 
% cellF2 = normr(cellF2); 

% 
% save(strcat('smtfF1F2.mat'),'cellF1F2','holdName')
% save(strcat('smtfF1.mat'),'cellF1','holdName')
% save(strcat('smtfF2.mat'),'cellF2','holdName')
% cd(typesFolder)
% save('sMTF Spot Size','sMTFSpotsOut')
%% choose cell types for area summ curve
% cellTypeLogical = contains(holdName,["ON Parasol","OFFParasol","ONSmooth","OFFSmooth"]);
cellTypeLogical = contains(holdName,["ON Parasol"]);
plot(areaRadii(cellTypeLogical,1:18)',areaCurve(cellTypeLogical,1:18)')
legend(holdName(cellTypeLogical))


%%
% 
% tsneF1F2=tsne(cellF1F2);
cellF1(isnan(cellF1))=0;
tsneF1=tsne(cellF1);
% tsneF2=tsne(cellF2);
figure
% gscatter(tsneF1F2(~strcmp(holdName,'Unknown'),1),tsneF1F2(~strcmp(holdName,'Unknown'),2),holdName(~strcmp(holdName,'Unknown')),tester)
% figure
gscatter(tsneF1(~strcmp(holdName,'Unknown'),1),tsneF1(~strcmp(holdName,'Unknown'),2),holdName(~strcmp(holdName,'Unknown')),colors2)
% figure
% gscatter(tsneF2(~strcmp(holdName,'Unknown'),1),tsneF2(~strcmp(holdName,'Unknown'),2),holdName(~strcmp(holdName,'Unknown')),tester)
% figure
% gscatter(tsneF1F2(:,1),tsneF1F2(:,2),holdName,tester)
% figure
% gscatter(tsneF1(:,1),tsneF1(:,2),holdName,tester)
% figure
% gscatter(tsneF2(:,1),tsneF2(:,2),holdName,tester)


%% Kmeans sort
celltoAnalyze = ["OFFParasol","OFFSmooth"];

cellnameList2 = cellnameList;

kmeansGroups2 = kmeansGroups;
kmeansGroups2(contains(cellnameList2,celltoAnalyze)==0,:)=[];
cellnameList2(contains(cellnameList2,celltoAnalyze)==0)=[];
cellnameList2(sum(kmeansGroups2(:,1:3),2)<5)=[];
kmeansGroups2(sum(kmeansGroups2(:,1:3),2)<5,:)=[]
% kmeansGroups2(sum(kmeansGroups2(strcmp(cellnameList,'ONParasol')==0,1:3),2)>5,:)=[];

%% mistake eraser

%     cd(strcat(exportFolder,fileIndex(:,1:4,d),'_',fileIndex(:,5:8,d)))
    splitFactors = ["noiseClass","epochGroupLabel","frameDwell","backgroundClass"];
    load(strcat(fileIndex(:,1:11,d),'.mat'))
    spikeEpochs=epochs;
    load(strcat(fileIndex(:,1:11,d),'_FT.mat'))
    frameTs = epochs;
    saveLabel = 'smTF';
    [splitCell,indexHolder,spikingData,frameTimings,metaData] = makeData(spikeEpochs,frameTs,protocolToAnalyze,splitFactors,6,recordingType);
    
%     cd(typeDir(1).folder)
    save(strcat(fileIndex(:,1:11,d),'_',saveLabel,'.mat'),'splitCell','indexHolder','frameTimings','-append')
