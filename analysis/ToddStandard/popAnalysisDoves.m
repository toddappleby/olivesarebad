%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'Doves Movie';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[]; %removes macOS hidden stuff
%%
 timings = [1250,6000,500]; %same for these
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
f1collection=[];
cellF1F2=[];
cellF1=[];
cellF2=[];
imageMeanAcrossCells=[];
pickParams = ["0" "2000" "50"];
sumsort = [];



% size(folderSet,3)
for c = 1:size(folderSet,3)
    
    cd(strcat(typesFolder,folderSet(:,:,c)))
    
    cellType = strtrim(folderSet(:,:,c))
    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    cellTypeData = [];
    
holdY=[];
        for d = 1:size(fileIndex,3)
      

         load(strtrim(fileIndex(:,:,d)))
         %Where params are most appropriate?  prob the ones I ran the most:
         sumsort=[];
         clear labelHolder
 for s = 1:size(indexHolder,2)
     labelHolder(s,:) = string(indexHolder{1,s});
 end
 
     for g = 1:size(labelHolder,1) % contains doesn't quite work when I have inner mask (mixes full field w/ surround conditions) so have this clumsy workaround
       sumsort(g,:) = sum(strcmp(labelHolder(g,:),pickParams),2);
     end
       if sumsort~=3
           continue
       end
       
       indicesWanted = indexHolder{2,sumsort==3};
        
       for p = 1:length(indicesWanted)
           psthImage(p,:) = psth(spikingData(indicesWanted(p),:),10,sampleRate,1);
       end
      
       averageForImage = mean(psthImage,1);
%        averageForImage = averageForImage/max(averageForImage);
       
       if isnan(averageForImage)
            continue
        end
       
       imageMeanAcrossCells = [imageMeanAcrossCells; averageForImage]; 

        
        holdName= [string(holdName); string(cellType)];

            
            
     

          
   
  
        end
        
        
       
        
        
        
        
end
% 
 cd(typesFolder)
 cd ..

% tester = imageMeanAcrossCells(strcmp(holdName,"Broad Thorny"),:);

save(strcat("dovesmovie_",pickParams(3)),'imageMeanAcrossCells')

%%
% 
load('/Users/reals/Documents/PhD 2021/ClarinetExports/Doves Movie/colors.mat')
figure
% tsneF1F2=tsne(cellF1F2);
% cellF1(isnan(cellF1))=0;
clear tsne1
tsne1=tsne(imageMeanAcrossCells);
% tsne1 = tsne(holder);
% tsneF2=tsne(cellF2);



figure
gscatter(tsne1(~strcmp(holdName,'Unknown'),1),tsne1(~strcmp(holdName,'Unknown'),2),holdName(~strcmp(holdName,'Unknown')),colors2)
figure
gscatter(tsne1(:,1),tsne1(:,2),holdName,colors)
figure
gscatter(tsne1(strcmp(holdName,'ON Parasol'),1),tsne1(strcmp(holdName,'ON Parasol'),2),holdName(strcmp(holdName,'ON Parasol')),colors2)
% figure
% gscatter(tsneF1(~strcmp(holdName,'Unknown'),1),tsneF1(~strcmp(holdName,'Unknown'),2),holdName(~strcmp(holdName,'Unknown')),tester)
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
    saveLabel = 'MN';
    [splitCell,indexHolder,spikingData,frameTimings,metaData] = makeData(spikeEpochs,frameTs,protocolToAnalyze,splitFactors,6,recordingType);
    
%     cd(typeDir(1).folder)
    save(strcat(fileIndex(:,1:11,d),'_','MN','.mat'),'splitCell','indexHolder','frameTimings','-append')
