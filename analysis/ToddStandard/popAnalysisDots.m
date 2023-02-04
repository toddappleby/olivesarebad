%super simple OMSDots .. 


%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'Object Motion Dots';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%celltypes must only contain type folders,  
%create array of folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[]; %removes hidden stuff
%% 
 typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    meanSC=[];
    semSC =[];
    for j = 3:size(fileSet,3)
        load(strcat(exportFolder,protocolToAnalyze,'/','celltypes/','ON OS/',fileSet(:,:,j)))
        
        uniqueSC = unique(splitCell{2,4});
        
        for h = 1:length(uniqueSC)
            
            scLogical = double(splitCell{2,4})==double(uniqueSC(h));
            meanSC(j-2,h) = mean(sum(spikingData(scLogical,:),2));
            semSC(j-2,h) = sem(sum(spikingData(scLogical,:),2));
            
            
            
        end
            meanSC(j-2,:) = meanSC(j-2,:)/max(meanSC(j-2,:));
        
    end
    
    
    

    
    plot(mean(meanSC))
    hold on
    plot(sem(meanSC))
    