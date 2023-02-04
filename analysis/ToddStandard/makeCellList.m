% get cell names and types from protocol-separated files

% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'S MT Fspot';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
% folderSet(:,:,1:3)=[]; %removes macOS hidden stuff
folderSet(:,:,1:2)=[];    
SMTFList(1,2)=string();
for a  = 1:size(folderSet,3)
    
    cd(strcat(typesFolder,folderSet(:,:,a))) 
    
    cellType = strtrim(folderSet(:,:,a));
    
    typeDir = dir;
    fileSet = char(string({typeDir.name})); 
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
       for b = 1:size(fileIndex,3)
       
          cellName = strtrim(fileIndex(:,1:11,b));
          
          SMTFList= [SMTFList; string(cellName) string(cellType)];
       
       end
        
end
SMTFList(1,:)=[];
cd(strcat(exportFolder,protocolToAnalyze))
save('SMTFCellList','SMTFList')

%% compare lists to find matching cells
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
compare1="S MT Fspot";
compare2="Doves Movie";
typesFolder = strcat(exportFolder,compare1);
cd(typesFolder)
load('SMTFCellList.mat')
typesFolder2 = strcat(exportFolder,compare2);
cd(typesFolder2)
load('DovesCellList.mat')

[lister,DovesIndex,SMTFIndex]=intersect(DovesList(:,1),SMTFList(:,1));

DovesMatch = DovesList(DovesIndex,:);
% SMTFMatch = SMTFList(SMTFIndex,:);
SMTFMatch = SMTFList(SMTFIndex,:);

DovesMatchSorted = sortrows(DovesMatch,2)
SMTFMatchSorted = sortrows(SMTFMatch,2)
cd(exportFolder)
save('SMTFandDovesMatch','DovesMatchSorted','SMTFMatchSorted')

%% find smtf with pulse not sinewave 

exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
cd(exportFolder)
load('SMTFandDovesMatch.mat')
cd(strcat([exportFolder,'/S MT Fspot/celltypes']))

dataDir = dir;
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[];    

% SMTFMatchSorted = sortrows(SMTFMatchSorted,2);
uniqueSMTF=unique(SMTFMatchSorted(:,2));
currentProtocol = "_smTF.mat";

cellHolder=cell(2,4);
count = 0; 

%%
for c = 1:length(unique(SMTFMatchSorted(:,2)))
    
    typeToFind = uniqueSMTF(c);
    
    cellList = SMTFMatchSorted(strcmp(typeToFind,SMTFMatchSorted(:,2)),1);
    
    cd(typeToFind)
    
    for d = 1:length(cellList)
        count=count+1
        load(strcat(cellList(d),currentProtocol))
         avgF1=[];
         avgF2=[];
         if sum(contains(splitCell{2,3},"sinewave"))>0
             if sum(contains(splitCell{2,2},"spot"))>0
                 matrixSplit = reshape([indexHolder{1,:}],[3,size(indexHolder,2)])';
                 logicalMatrix = contains(matrixSplit,["spot" "sinewave"]);
                 outI = find(sum(logicalMatrix,2)==2);
             else
                  count = count -1;
                 continue
                
             end
         else
             count = count-1;
             continue
             
         end
                 
         wantedIndices = indexHolder{2,outI};
         Radii = splitCell{2,4}(wantedIndices);
         Radii = str2double(Radii);
         
         uRadii = unique(Radii);
      for toFindMultiples = 1:length(unique(Radii))
          
          indexofIndex = find(Radii==uRadii(toFindMultiples));
          finalIndices=wantedIndices(indexofIndex);
            collectF1=[];
            collectF2=[];
            collectPhase=[];
            
         for multiplesCounter = 1:length(finalIndices)
                binnedData = BinSpikeRate(spikingData(finalIndices(multiplesCounter),5001:30000),100,10000)';
                [F, phase] = frequencyModulation(binnedData, ...
                100, 2, 'avg', 1:2, []);
                
            collectF1 = [collectF1; F(1)];
            collectF2 = [collectF2; F(2)];
            collectPhase = [collectPhase; phase];
         end
         
         avgF1(toFindMultiples) = mean(collectF1);
         avgF2(toFindMultiples) = mean(collectF2);
         phaser(toFindMultiples,:) = mean(collectPhase,1);
         
      end

      cellHolder{count,1}=typeToFind;
      cellHolder{count,2}=cellList(d);
      cellHolder{count,3}=avgF1;
      cellHolder{count,4}=avgF2;
      cellHolder{count,5}=uRadii;
         
    end

    cd .. 
    
end

%% 
figure(98)
counter=0;
for z = 20
    counter=counter+1;
    plot(cellHolder{z,5},cellHolder{z,3})
    hold on
    forF1Mean(counter,1:18) = cellHolder{z,3} 
    z
    
end
figure
plot(uRadii,mean(forF1Mean,1))
% hold on
% plot(sem(forF1Mean))
makeAxisStruct(gca,strtrim(['smtfspotMean' '_Recursive']))
