%% 
exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
protocolToAnalyze = 'Motion and Noise';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:3)=[]; %removes macOS hidden stuff
%%
kmeansGroups=[];
cellnameList = [];
errorGroups = [];
for c = 1:size(folderSet,3)
    
    cd(strcat(typesFolder,folderSet(:,:,c)))
    
    cellType = strtrim(folderSet(:,:,c));
    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    cellTypeData = [];
    
    
        for d = 1:size(fileIndex,3)

         load(strtrim(fileIndex(:,:,d)))
         
          cd(strcat(exportFolder,fileIndex(:,1:4,d),'_',fileIndex(:,5:8,d)))
    splitFactors = ["noiseClass","epochGroupLabel","frameDwell","backgroundClass"];
    load(strcat(fileIndex(:,1:11,d),'.mat'))
    saveLabel = 'MN';
    [splitCell,indexHolder,spikingData,metaData] = makeData(epochs,protocolToAnalyze,splitFactors,6,recordingType);
    
    cd(typeDir(1).folder)
    save(strcat(fileIndex(:,1:11,d),'_','MN','.mat'),'splitCell','indexHolder','-append')

             % "contrast","centerBarOrientation","surroundBarOrientation"
                % "centerBarWidth","numberOfCenterBars","backgroundClass","centerClass,"onlineAnalysis""
%         combosWanted = ["1","90","90","sequential","sequential"
%                        "1","90","90","sequential","random"
%                        "1","90","90","sequential","sequential180"
%                        "1","90","90","random","sequential"
%                         "1","90","90","random","random"
%                         "1","90","90","random","sequential180"]; %yes, this is a little redundant w/ index holder
% 


                p = [indexHolder{1,:}];    
                cellToAnalyze=[];
                   
                    for f = 1:size(combosWanted,1) %rows

                        
                        pLogic = ismember(p,combosWanted(f,:));
                        pLogic = reshape(pLogic,size(splitCell,2),size(indexHolder,2)); %reshape lists column wise (i.e. 1 4 7; 2 5 8; 3 6 9)
                        pLogic = pLogic';
                        contractToMatch = logical(sum(pLogic,1));
                        pMatrix = reshape(p,size(splitCell,2),size(indexHolder,2));
                        pMatrix = pMatrix'; 
                        pMatrix = pMatrix(:,contractToMatch); %not sure this is the greatest coding moment in history .. 

                            if size(pMatrix,2) ~= size(combosWanted,2)
                                cellToAnalyze = [];
                                 continue
                            end
                            
                        pIndex = ismember(pMatrix,combosWanted(f,:),'rows');
                            

                        cellTransfer = indexHolder(2,pIndex)';
                        finalInd = vertcat(cellTransfer{:});
                        finalInd = finalInd';


                        cellToAnalyze{f,1} = join(combosWanted(f,:));
                        cellToAnalyze{f,2} = finalInd;

                    end
     
            if isempty(cellToAnalyze)
                continue
            end
                    
            cellToAnalyze = cellToAnalyze';
            
            
            for g = 1:size(cellToAnalyze,2)
            cellTypeData(d,g) = mean(sum(spikingData(cellToAnalyze{2,g},11800:15000),2));
            cellTypeDataError(d,g) = sem(sum(spikingData(cellToAnalyze{2,g},11800:15000),2));
            end
            
  
      cellnameList = [cellnameList; string(cellType)];    
      

    end
    cellTypeData(isnan(cellTypeData))=0;
    cellTypeDataError(isnan(cellTypeDataError))=0;
    
    
           
    kmeansGroups = [kmeansGroups; cellTypeData];
    errorGroups = [errorGroups; cellTypeDataError];
%     kmeansGroups = kmeansGroups ./ max(abs(kmeansGroups),[],2);
    
   kmeansGroups = kmeansGroups(sum(kmeansGroups,2)>0,:);
   errorGroups = errorGroups(sum(kmeansGroups,2)>0,:);
   cellnameList = cellnameList(sum(kmeansGroups,2)>0,:);
    
    cellTypeData = cellTypeData(sum(cellTypeData,2)>0,:); % these two lines do the same thing 
 %get rid of cells that weren't run on selected parameters
    
    dataCell{1,c}=cellType;
    dataCell{2,c}=mean(cellTypeData,1);
    
    errorCell{1,c} = cellType;
    errorCell{2,c} = mean(cellTypeDataError,1);
    
    
end

%% graphing
figure(3)
for h = 1:length(dataCell)
    
    subplot(3,3,h)
    hold on
    bar(dataCell{2,h})
    ylabel('spike count')
%     set(gca,'xticklabel',categorical({'SS','SR','S180','RS','RR','R180'}))
    title(dataCell{1,h})
    er = errorbar(dataCell{2,h},errorCell{2,h});
    er.Color = [0 0 0];
    er.LineStyle = 'none';

end
sgtitle('Summary: Cell Type Responses to Dark Bars')
legend('Left to right (center/surround): motion/motion,random/motion,180/motion,motion/non-motion,random/non-motion,180/non-motion')
%% Kmeans sort
celltoAnalyze = ["OFFParasol","OFFSmooth"];

cellnameList2 = cellnameList;

kmeansGroups2 = kmeansGroups;
kmeansGroups2(contains(cellnameList2,celltoAnalyze)==0,:)=[];
cellnameList2(contains(cellnameList2,celltoAnalyze)==0)=[];
cellnameList2(sum(kmeansGroups2(:,1:3),2)<5)=[];
kmeansGroups2(sum(kmeansGroups2(:,1:3),2)<5,:)=[]
% kmeansGroups2(sum(kmeansGroups2(strcmp(cellnameList,'ONParasol')==0,1:3),2)>5,:)=[];
%%
cell1 = "OFFParasol"
cell2 = "OFFSmooth"

normData = normr(kmeansGroups2(:,1:3));
normDataRsurr = normr(kmeansGroups2(:,4:6));

% [idx,centroids] = kmeans(kmeansGroups(:,1:3),6);
[idx1,centroids] = kmeans(normData,3);



figure(2);scatter3(normData(strcmp(cellnameList2,cell1),1), normData(strcmp(cellnameList2,cell1),2), normData(strcmp(cellnameList2,cell1),3),40,categorical(cellnameList2(strcmp(cellnameList2,cell1))), 'filled','g')
hold on;figure(2);scatter3(normData(strcmp(cellnameList2,cell2),1), normData(strcmp(cellnameList2,cell2),2), normData(strcmp(cellnameList2,cell2),3),40,categorical(cellnameList2(strcmp(cellnameList2,cell2))), 'filled','b')

% figure(2);scatter3(normData(:,1), normData(:,2), normData(:,3), 40,idx1, 'filled')

% xlabel('seq');
% ylabel('rand');
% zlabel('180');
xlim([0 1])
ylim([0 1])
zlim([0 1])

% text(normData(:,1), normData(:,2), normData(:,3),cellnameList2)

[idx2,centroids] = kmeans(normDataRsurr,3);
figure(2);scatter3(normDataRsurr(strcmp(cellnameList2,cell1),1), normDataRsurr(strcmp(cellnameList2,cell1),2), normDataRsurr(strcmp(cellnameList2,cell1),3),40,categorical(cellnameList2(strcmp(cellnameList2,cell1))), '+','g')

figure(2);scatter3(normDataRsurr(strcmp(cellnameList2,cell2),1), normDataRsurr(strcmp(cellnameList2,cell2),2), normDataRsurr(strcmp(cellnameList2,cell2),3),40,categorical(cellnameList2(strcmp(cellnameList2,cell2))), '+','b')

% figure(2);hold on;scatter3(normDataRsurr(:,1),normDataRsurr(:,2),normDataRsurr(:,3),40,categorical(cellnameList2),'+')
legend(strcat(cell1,' (Motion Surround)'),strcat(cell2,' (Motion Surround)'),strcat(cell1,' (Random Surround)'),strcat(cell2,' (Random Surround)'))


xlabel('sequential');

ylabel('random');
zlabel('180');
xlim([0 1])
ylim([0 1])
zlim([0 1])
% colormap(gca,'turbo')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9, 9], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
title('OFF Parasol & OFF Smooth Cell Responses for Seq vs Random vs 180 center conditions w/ each surround')
% text(normDataRsurr(:,1), normDataRsurr(:,2), normDataRsurr(:,3),cellnameList2)
%% thing
cd('/Users/toddappleby/Documents/Data/Clarinet Exports/Motion Center Surround')
cellnameList2 = cellnameList2';
save('MCSKmeans','normData','normDataRsurr','idx1','idx2')
save('MCSlabels','cellnameList2')
%% mistake eraser

 cd(strcat(exportFolder,fileIndex(:,1:4,d),'_',fileIndex(:,5:8,d)))
    splitFactors = ["noiseClass","epochGroupLabel","frameDwell","backgroundClass"];
    load(strcat(fileIndex(:,1:11,d),'.mat'))
    saveLabel = 'MN';
    [splitCell,indexHolder,spikingData,metaData] = makeData(epochs,protocolToAnalyze,splitFactors,6,recordingType);
    
    cd(typeDir(1).folder)
    save(strcat(fileIndex(:,1:11,d),'_','MCS','.mat'),'splitCell','indexHolder','-append')

