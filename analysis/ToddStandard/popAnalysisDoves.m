%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'Doves Movie';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%celltypes must only contain type folders,  
%create array of folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[]; %removes hidden stuff
%%
 timings = [1250,6000,500]; %same for these
 sampleRate = 10000;

stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;

kmeansGroups=[];
cellnameList = [];
errorGroups = [];
recordingType = 'extracellular';
binRate = 1e3;
holdName=[];
imageMeanAcrossCells=[];
pickParams = ["0" "2000" "30"];
sumsort = [];
holdFile=[];


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
           psthImage(p,:) = binSpikeCount(spikingData(indicesWanted(p),:),60,sampleRate);
%            psthImage(p,:) = psth(spikingData(indicesWanted(p),:),10,sampleRate,1);
       end
      
       averageForImage = mean(psthImage,1);
%        averageForImage = averageForImage/max(averageForImage);
       
       if isnan(averageForImage)
            continue
       end
        
       
       imageMeanAcrossCells = [imageMeanAcrossCells; averageForImage]; 

        
        holdName= [string(holdName); string(cellType)];
        holdFile = [string(holdFile); string(fileIndex(:,1:11,d))];
        end
end
% 
 cd(typesFolder)
 cd ..


% tester = imageMeanAcrossCells(strcmp(holdName,"Broad Thorny"),:);

% save(strcat("dovesmovie_",pickParams(3)),'imageMeanAcrossCells')
%% SVD by cell type
nameList = unique(holdName);
svdApprox = [];
for eachType = 1:size(unique(holdName),1)
    typeDataSVD = imageMeanAcrossCells(strcmp(holdName,nameList(eachType)),:);
    coMatrix = typeDataSVD' * typeDataSVD;
    [U,S,V] = svd(coMatrix);
    rank = 2;
    recreate = U(:,rank)*S(rank,rank)*V(:,rank)';
    
    svdApprox = [svdApprox; typeDataSVD*recreate;]; %#ok<AGROW>
    
end
%% eig it
nameList = unique(holdName);
eigApprox = [];
for eachType = 1:size(unique(holdName),1)
    
    typeDataSVD = imageMeanAcrossCells(strcmp(holdName,nameList(eachType)),:);
    [m,n]=size(typeDataSVD); 
    mn=mean(typeDataSVD,2);
%     typeDataSVD=typeDataSVD-repmat(mn,1,n);
%     coMatrix =(1/(n-1))*typeDataSVD' * typeDataSVD;
coMatrix =typeDataSVD' * typeDataSVD;
    [V,D] = eig(coMatrix);
    lambda=diag(D); 
    [sortthing,m_arrange]=sort(-1*lambda);  %eigenvalues come out backwards 
    V=V(:,m_arrange);  
    tester=V'*coMatrix;
    
    eigApprox = [eigApprox; tester(1,:).*typeDataSVD;]; %#ok<AGROW>
    
end
%% lol
imageMeanAcrossCells(92,:)=r(1,:)
holdName(92)="RecursiveMike"
%% tsne
% 
 load('/Users/reals/Documents/PhD 2021/ClarinetExports/Doves Movie/colors.mat')

clear tsne1
tsne1=tsne(imageMeanAcrossCells);




figure
gscatter(tsne1(~strcmp(holdName,'Unknown'),1),tsne1(~strcmp(holdName,'Unknown'),2),holdName(~strcmp(holdName,'Unknown')),colors2)
figure
gscatter(tsne1(:,1),tsne1(:,2),holdName,colors)
figure
gscatter(tsne1(strcmp(holdName,'ON Parasol'),1),tsne1(strcmp(holdName,'ON Parasol'),2),holdName(strcmp(holdName,'ON Parasol')),colors2)
%% look at responses to single cell type
typeWanted = "OFF Parasol";
typeData = imageMeanAcrossCells(strcmp(holdName,typeWanted),:);

% typeWanted2 = "OFF Parasol";
% typeData2 = imageMeanAcrossCells(strcmp(holdName,typeWanted2),:);
figure
for cellNum = 1:size(typeData,1)
    subplot(size(typeData,1),1,cellNum)
  
    plot(typeData(cellNum,:))
%     hold on
%     plot(typeData2(cellNum+2,:))
end
% legend(typeWanted,typeWanted2)
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
