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

protocolToFind = "Doves Movie";


%%
cellList = [];
for f = 1:size(folderIndex,3)
    
    
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
mkdir(protocolToFind)
cd(strcat('/Users/reals/Documents/PhD 2021/ClarinetExports/',protocolToFind))
save(strtrim(protocolToFind),'cellList')
%% and now we parse spikes 
protocolToOpen = 'Doves Movie';
cd(strcat(exportFolder,protocolToOpen))
load(strcat(protocolToOpen,'.mat'));
recordingType = 'extracellular';
splitFactors = ["maskDiameter","apertureDiameter","stimulusIndex"];
saveLabel = 'Doves';
cellList=cellList';


%% ok, NOW we parse spikes

    for b = 99:length(cellList)
    cellName = char(cellList(b));
    
    %date (year/month/day) by naming convention 
    cellYear = cellName(1:4);
    cellDate = cellName(5:8);
    load(strcat(exportFolder,cellYear,'_',cellDate,'/',cellName))
    
    [splitCell,indexHolder,spikingData,~,metaData,~] = makeData(epochs,[],protocolToOpen,splitFactors,3,recordingType);
    
save(strcat(cellName(1:11),'_','MN'),'splitCell','indexHolder','spikingData','metaData')
% save(strcat(cellName(1:11),'_','MCS'),'splitCell','indexHolder','metaData','-append') %index fix haha oops (no spikes)
end

    
%147, use '-append' to change indexholder