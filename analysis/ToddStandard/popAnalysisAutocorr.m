%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/';
protocolToAnalyze = 'Motion And Noise';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:3)=[]; %removes macOS hidden stuff
%%
 timings = [250,10000,250]; %same for these
 nonlinearityBins = 100;
kmeansGroups=[];
cellnameList = [];
errorGroups = [];
recordingType = 'extracellular';
binRate = 1e3;
yCrossCells=[];
holdName=[];
% size(folderSet,3)
for c = 1:length(folderSet)
    
    cd(strcat(typesFolder,folderSet(:,:,c)))
    
    cellType = strtrim(folderSet(:,:,c))
    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    cellTypeData = [];
    
holdY=[];
        for d = 1:size(fileIndex,3)
           outI = []; 

         load(strtrim(fileIndex(:,:,d)))
         %Where params are most appropriate?  prob the ones I ran the most:
%          clear indLength
%             for e = 1: length(indexHolder)
%                 indLength(e) = length(indexHolder{2,e});
%             end
% btwn is new

 sumsortSeq=[];
  sumsortRand=[];
   sumsortStatic=[];
         clear labelHolder
 for s = 1:size(indexHolder,2)
     labelHolder(s,:) = string(indexHolder{1,s});
 end

       sumsortSeq = sum(contains(labelHolder,["binary" "1" "sequential"]),2);
       sumsortRand = sum(contains(labelHolder,["binary" "1" "random"]),2);
       sumsortStatic = sum(contains(labelHolder,["binary" "1" "stationary"]),2);
       if sumsortStatic~=3
           continue
       end
       
%        indicesWanted = indexHolder{2,sumsort==2};
%  indicesWanted = indexHolder{2,sumsort==3};
%btwn is new
            
%             maxI = maxk(indLength,3); %top 3 conditions -- must be best options for seq rand & static
%             logI = ismember(string(indLength),string(maxI));
%             outI = find(logI==1); %index holder indices of top options
            
            %which indices correspond to which condition?
%             clear holdLabel
%             for ee = 1:length(outI)
%                 holdLabel(ee) = string(indexHolder{1,outI(ee)}(4));
%             end
            
%             seqIndex = indexHolder{2,outI(strcmp(holdLabel,"sequential"))};
%             randIndex = indexHolder{2,outI(strcmp(holdLabel,"random"))};
%             staticIndex = indexHolder{2,outI(strcmp(holdLabel,"stationary"))};

            seqIndex = indexHolder{2,sumsortSeq==3};
            randIndex = indexHolder{2,sumsortRand==3};
            staticIndex = indexHolder{2,sumsortStatic==3};


            
            staticSpikes = spikingData(staticIndex,:);
            
            G= find(staticSpikes(:) >0);
            isi = diff(G);
            isi=isi/10;
            [y x] = getSpikeAutocorrelation(isi);
%             y=y/max(y);
            

          
       
            
%             
%             prePts = 250 * 1e-3 * binRate;
%             stimPts = 10000 * 1e-3 * binRate;
%             tailPts = 250 * 1e-3 * binRate; 
            
            holdY(d,:)=y';
            yCrossCells = [yCrossCells; y];
             holdName = [holdName string(cellType)];
        end
        
        
       
        
        
        figure

plot(mean(holdY,1))
        title(cellType)
       
saveY=[];
saveY=mean(holdY,1)';
save(strcat('AutoCorr','_',cellType,'.mat'),'saveY')
        
end

cd(typesFolder)
save(strcat('AutoCorrfortSNE.mat'),'yCrossCells','holdName')

%tsne
%% tsne
[i,j]=find(~isnan(yCrossCells))
sorter=unique(i);
holdName=holdName(sorter);
yCrossCells=yCrossCells(sorter,:);
ACtSNE=tsne(yCrossCells);
gscatter(ACtSNE(:,1),ACtSNE(:,2),holdName')

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
