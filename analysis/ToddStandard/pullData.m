%% 
%This script should allow opening of .mat files, looking at protocols, grabbing 
%spikes from desired protocols, then saving as new .mats.  Should work in
%any data directory w/ consistent naming scheme.  Should work for data
%aggregation with change to search for protocol in filename.  filename
%probably should contain protocol and celltype along with date & cell num
%% setup
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/'
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/'
cd(exportFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderIndex = folderSet(:,1:9,folderSet(:,1,:) == '2');% logical index in 3rd dim, 1:9 deletes spaces that formed for some reason at end of folder name after converting to char array

protocolToFind = "Motion And Noise";


%%
cellList = [];
for f = 73:size(folderIndex,3)
    
    
    cd(strcat(exportFolder,folderIndex(:,:,f)))
    cellDir = dir;
    cellSet = char(string({cellDir.name}));
    cellNames = cellSet(:,:,cellSet(:,1,:) == '2');
    
    folderIndex(:,:,f)    
    
        for a = 1:size(cellNames,3)
            if ~contains(strcat(strtrim(cellNames(:,:,a))),'_FT') %so doesn't process frame timing files (redundant)
                
                load(strcat(strtrim(cellNames(:,:,a))))

                
                uniqueProtocols = getProtocols(epochs);
                
                    if sum(strcmp(uniqueProtocols,protocolToFind))>0
                        cellList = [cellList string(strtrim(cellNames(:,:,a)))];
                    end
            end
        end
end
cd(exportFolder)


% below adds new directory if never pulled protocol, otherwise checks old
% list against new list for new files (note: only files after end file in old list to
% avoid processing abandoned [bad] cells again)
if ~exist(protocolToFind,'dir')
mkdir(protocolToFind)
else
cellList2=cellList;
load(strcat(protocolToFind,'.mat'))
cellListFinal = cellList2(find(strcmp(cellList2,cellList(length(cellList))))+1:end); %i don't know another way to find index for the 1 in the logical....  must be more clever way to do this 
cellList = [cellList cellListFinal];
end




cd(strcat('/Users/reals/Documents/PhD 2021/ClarinetExports/',protocolToFind))
save(strcat(strtrim(protocolToFind),'.mat'),'cellList')
%% and now we parse spikes // note: just load protocol mat file if doing corrections (don't run above section -- you can, but it's a time waste)
parseNewData = 1; %this is set up BACKWARDS for code segment after this! fuck!
protocolToOpen = 'Motion And Noise';
cd(strcat(exportFolder,protocolToOpen))
load(strcat(protocolToOpen,'.mat'));
recordingType = 'extracellular';
% splitFactors = ["radius","splitContrasts","contrast","spaceConstant"];
% splitFactors = ["maskDiameter","apertureDiameter","stimulusIndex"]; %DOVES
splitFactors = ["noiseClass","epochGroupLabel","frameDwell","apertureRadius","barOrientation","backgroundClass"]; %MotionAndNoise 
% splitFactors = ["spotIntensity","backgroundIntensity","currentSpotSize"]; %Exp Spots I guess
% splitFactors = ["contrast","stimulusClass","temporalClass","radius"]; %smtf
% splitFactors = ["stimulusClass"]; %OMS
% splitFactors = ["intensity","backgroundIntensity","barSize","orientation"];
% 
saveLabel = 'MN';
if parseNewData
    cellList = cellListFinal'; %parse only new cells 
else
    cellList=cellList'; %parse whole thang
end



%% ok, NOW we parse spikes
%
for b = 45:length(cellList)
    cellName = char(cellList(b));
    
    %date (year/month/day) by naming convention 
    cellYear = cellName(1:4);
    cellDate = cellName(5:8);
    cellName = strcat(cellName(1:11),'_FT','.mat');
    load(strcat(exportFolder,cellYear,'_',cellDate,'/',cellName)) %for m and n 
    frameTs=epochs;
    load(strcat(exportFolder,cellYear,'_',cellDate,'/',char(cellList(b))))
    
    
[splitCell,indexHolder,spikingData,frameTimings,metaData,seed] = makeData(epochs,frameTs,protocolToOpen,splitFactors,5,recordingType,parseNewData);
    
%     [splitCell,indexHolder,spikingData,~,metaData,~] = makeData(epochs,[],protocolToOpen,splitFactors,6,recordingType,parseNewData);
    
save(strcat(cellName(1:11),'_',saveLabel),'splitCell','indexHolder','spikingData','metaData')
% save(strcat(cellName(1:11),'_',saveLabel),'splitCell','indexHolder','metaData')
% save(strcat(cellName(1:11),'_',saveLabel),'seed','frameTimings','splitCell','indexHolder','-append') %index fix haha oops (no spikes)
end

%use '-append' to change indexholder
 %% really annoying behavior fix (late) -- redoing splitCell and indexHolder for each cell type manually (have to navigate to folder then run this section followed by above section)
%  oldCellList=cellList; do this first then comment out
 thing = dir;
 fileToFix=string();
for i = 3:length(thing)-1
    fileToFix(i-2) = string(thing(i).name(1:11));
    
end
 fileToFix=fileToFix';
 cellList = oldCellList(contains(oldCellList,fileToFix));