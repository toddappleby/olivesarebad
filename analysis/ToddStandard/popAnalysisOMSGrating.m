%%OMS Texture 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
% exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
exportFolder = 'C:\Users\todda\Documents\Primate Data\ClarinetExports\';
protocolToAnalyze = 'Object Motion Grating';
typesFolder = strcat(exportFolder,protocolToAnalyze,'\','celltypes\');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[]; %removes macOS hidden stuff
%%
  %same for these
sampleRate = 10000; 
nonlinearityBins = 100;
kmeansGroups=[];
cellnameList = [];
errorGroups = []; 
recordingType = 'extracellular';
binRate = 1e3; 
yCrossCells=[];
holdType=[];
expSpotsOut=cell(0);
count=0;

centerMean = [];
diffMean = [];
globalMean = [];
surroundMean = [];


centerError = [];
diffError = [];
globalError = [];
surroundError = [];

% size(folderSet,3)
for c = 2:2            

cd(strcat(typesFolder,folderSet(:,:,c)))

cellType = strtrim(folderSet(:,:,c))

typeDir = dir;
fileSet = char(string({typeDir.name}));
fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');

centerMeanType = [];
diffMeanType = [];
globalMeanType = [];
surroundMeanType = [];

centerForShow=[];
diffForShow=[];
globalForShow=[];
surroundForShow=[];

        for d = 1:size(fileIndex,3)
     
        allCenterSpikes= [];
        allCenterSpikesPSTH = [];
        allDiffSpikes = [];
        allDiffSpikesPSTH = [];
        allGlobalSpikes = [];
        allGlobalSpikesPSTH = [];
        allDriftSpikes = [];
        allDriftSpikesPSTH = [];



         load(strtrim(fileIndex(:,:,d)))
         
         cellDate = fileIndex(:,1:11,d);
         
         % get timings for each cell just in case diff numbers used 
         timings = [metaData.preTime,metaData.stimTime,metaData.tailTime];
         stimStart = (timings(1)*1e-3)*sampleRate+1;
         stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
         stimOff = (timings(1)+timings(2)+timings(3))*10;
         % timings done
         
         %to match mike: differential, global, center, drift   
         % ~!!!!!! SWITCH INDEX HOLDER NUMBERS BETWEEN TEXTURE AND GRATING

         allCenterSpikes = spikingData(indexHolder{2,4},stimStart:stimEnd);

         for aa = 1:size(allCenterSpikes,1)
    
         allCenterSpikesPSTH(aa,:) = psth(allCenterSpikes(aa,:),6+2/3,10000,1);

         end

         centerPSTHmean = mean(allCenterSpikesPSTH,1);

         %diff or eyeobj
         allDiffSpikes = spikingData(indexHolder{2,2},stimStart:stimEnd);

         for aa = 1:size(allDiffSpikes,1)
    
         allDiffSpikesPSTH(aa,:) = psth(allDiffSpikes(aa,:),6+2/3,10000,1);

         end

         diffPSTHmean = mean(allDiffSpikesPSTH,1);

        %global or eye
         allGlobalSpikes = spikingData(indexHolder{2,1},stimStart:stimEnd);

         for aa = 1:size(allGlobalSpikes,1)
    
         allGlobalSpikesPSTH(aa,:) = psth(allGlobalSpikes(aa,:),6+2/3,10000,1);

         end

         globalPSTHmean = mean(allGlobalSpikesPSTH,1);

         %surround or drift 
         allDriftSpikes = spikingData(indexHolder{2,3},stimStart:stimEnd);

         for aa = 1:size(allDriftSpikes,1)
    
         allDriftSpikesPSTH(aa,:) = psth(allDriftSpikes(aa,:),6+2/3,10000,1);

         end

         driftPSTHmean = mean(allDriftSpikesPSTH,1);


avgData = [diffPSTHmean; globalPSTHmean; centerPSTHmean; driftPSTHmean];
classLabel = [indexHolder{1,2}; indexHolder{1,1}; indexHolder{1,4}; indexHolder{1,3}];



save(strcat('C:\Users\todda\Documents\Primate Data\Saves\',cellType,'_',cellDate,'_Todd.mat'),'avgData','classLabel')


%          centerSpikes = mean(sum(spikingData(indexHolder{2,1},stimStart:stimEnd),2));
%          diffSpikes = mean(sum(spikingData(indexHolder{2,2},stimStart:stimEnd),2));
%          globalSpikes = mean(sum(spikingData(indexHolder{2,3},stimStart:stimEnd),2));
%          surroundSpikes = mean(sum(spikingData(indexHolder{2,4},stimStart:stimEnd),2));
% 
% 
%          
%          for z = 1:length(indexHolder{2,3})
%          centerForShow(d,:) = binSpikeCount(spikingData(indexHolder{2,1}(z),stimStart:stimEnd),60,10000);
%          diffForShow(d,:)  = binSpikeCount(spikingData(indexHolder{2,2}(z),stimStart:stimEnd),60,10000);
%          globalForShow(d,:)  = binSpikeCount(spikingData(indexHolder{2,3}(z),stimStart:stimEnd),60,10000);
%          surroundForShow(d,:)  = binSpikeCount(spikingData(indexHolder{2,4}(z),stimStart:stimEnd),60,10000);
%          end
%         
%          centerMeanType = [centerMeanType centerSpikes];
%          diffMeanType = [diffMeanType diffSpikes];
%          globalMeanType = [globalMeanType globalSpikes];
%          surroundMeanType = [surroundMeanType surroundSpikes];
%         
%          if d==20000
%              subplot(3,1,1)
%              plot(mean(centerForShow))
%              makeAxisStruct(gca,strtrim(['OMSGrating' '_OffPCenter']))
%               subplot(3,1,2)
%              plot(mean(diffForShow))
%              makeAxisStruct(gca,strtrim(['OMSGrating' '_OffPDiff']))
%               subplot(3,1,3)  
%              plot(mean(globalForShow))
%              makeAxisStruct(gca,strtrim(['OMSGrating' '_OffPGlobal']))
% 
%          
%          end
%         end
%        
%         fullDist=[];
%         fullDist = [centerMeanType; diffMeanType; globalMeanType; surroundMeanType];
%         normDist = fullDist./max(fullDist);
        
%         centerMean = [centerMean; mean(normDist(1,:),2)];
%         diffMean = [diffMean;  mean(normDist(2,:),2)];
%         globalMean = [globalMean;  mean(normDist(3,:),2)];
%         surroundMean = [surroundMean;  mean(normDist(4,:),2)];
%         
%         centerError = [centerError; sem(normDist(1,:),2)];
%         diffError = [diffError; sem(normDist(2,:),2)];
%         globalError = [globalError; sem(normDist(3,:),2)];
%         surroundError = [surroundError; sem(normDist(4,:),2)];
%         
%         centerMean = [centerMean; mean(centerMeanType)];
%         diffMean = [diffMean; mean(diffMeanType)];
%         globalMean = [globalMean; mean(globalMeanType)];
%         surroundMean = [surroundMean; mean(surroundMeanType)];
%         
%         centerError = [centerError; sem(centerMeanType)];
%         diffError = [diffError; sem(diffMeanType)];
%         globalError = [globalError; sem(globalMeanType)];
%         surroundError = [surroundError; sem(surroundMeanType)];
        
        
%          holdType = [string(holdType); string(cellType)];
       
        
        
        
        
        end
end

 cd(typesFolder)
%% plot it

plot([centerMean(2) diffMean(2) globalMean(2) surroundMean(2)])
makeAxisStruct(gca,strtrim(['OMSGratingMeans' '_BT']))
figure
plot([centerMean(3) diffMean(3) globalMean(3) surroundMean(3)])
makeAxisStruct(gca,strtrim(['OMSGratingMeans' '_OffP']))
figure
plot([centerMean(4) diffMean(4) globalMean(4) surroundMean(4)])
makeAxisStruct(gca,strtrim(['OMSGratingMeans' '_OffNarrowThorny']))
figure
plot([centerMean(5) diffMean(5) globalMean(5) surroundMean(5)])
makeAxisStruct(gca,strtrim(['OMSGratingMeans' '_Recursive']))

