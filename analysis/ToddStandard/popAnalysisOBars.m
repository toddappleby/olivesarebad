%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'Oriented Bars';
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
cellnameList = [];
recordingType = 'extracellular';
% pickParams = ["-1" "0.5" "100  1000"];
% pickParams = ["-1" "0.5" "50  1000"];
pickParams = ["1" "0" "100  1000"];
% pickParams = ["1" "0.5" "100  1000"];
binRate = 1e3;
yCrossCells=[];
holdName=[];
holdFile=[];
count=0;
areaCurve=[];
areaRadii=[];

% cellCollect = cell(size(folderSet,3),2);

% 
for c = 1:size(folderSet,3)
    
    cd(strcat(typesFolder,folderSet(:,:,c)))
    
    cellType = strtrim(folderSet(:,:,c))
    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    cellTypeData = [];
    


for d = 1:size(fileIndex,3)
meanOffNull=[];
meanOnNull=[];
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

 if strcmp(fileIndex(:,:,d),'20220223Ac2_OBars.mat')
     pickParams = ["-1" "0.5" "150  1000"];
 end
 
  if strcmp(fileIndex(:,:,d),'20210907Ac2_OBars.mat')
     pickParams = ["-1" "0.5" "50  1000"];
 end
 
  for g = 1:size(labelHolder,1) % contains doesn't quite work when I have inner mask (mixes full field w/ surround conditions) so have this clumsy workaround
   sumsort(g,:) = sum(strcmp(labelHolder(g,:),pickParams),2);
  end

   if sumsort~=3
       
       d
       continue
   end

    indicesWanted = indexHolder{2,sumsort==3};
   compileOrientations = [];
   offResp = [];
   onResp = [];
    for p = 1:length(indicesWanted)
        if double(pickParams(1)) < 0 
               offResp(p) = sum(spikingData(indicesWanted(p),5000:10000),2);
               onResp(p) = sum(spikingData(indicesWanted(p),10000:15000),2);
               compileOrientations(p) = splitCell{2,4}(indicesWanted(p)); 
        else 
               onResp(p) = sum(spikingData(indicesWanted(p),5000:10000),2);
               offResp(p) = sum(spikingData(indicesWanted(p),10000:15000),2);
               compileOrientations(p) = splitCell{2,4}(indicesWanted(p)); 
        end
    end
    uniques = unique(compileOrientations);
    meanOff=[];
    meanOn=[];
    for n = 1:length(uniques)
        meanOff(n) = mean(offResp(compileOrientations == uniques(n)));
        meanOn(n) = mean(onResp(compileOrientations == uniques(n)));
    end
          meanOffMax = max(meanOff);
        
        meanOffMaxOri = uniques(meanOff==meanOffMax);
        
        if meanOffMaxOri >= 90
            meanOffNullOri = meanOffMaxOri(1)-90;
        else
            meanOffNullOri = meanOffMaxOri(1)+90;
        end
        
        for q = -1:1 %find lowest val near null
         
            transferNull = meanOff(uniques==meanOffNullOri + (q*15));                  
          
            if ~isempty(transferNull)
                meanOffNull = [meanOffNull transferNull];
            end
        end
        meanOffNull=min(meanOffNull);
        %on
        meanOnMax = max(meanOn);
        
        meanOnMaxOri = uniques(meanOn==meanOnMax);
        
        if meanOnMaxOri >= 90
            meanOnNullOri = meanOnMaxOri(1)-90;
        else
            meanOnNullOri = meanOnMaxOri(1)+90;
        end
        
        for q = -1:1 %find lowest val near null
         
            transferNull = meanOn(uniques==meanOnNullOri + (q*15));                  
          
            if ~isempty(transferNull)
                meanOnNull = [meanOnNull transferNull];
            end
        end
          meanOnNull=min(meanOnNull);
        
% 
        
        
      
        
        
         cellCollect3{1,count} = string(cellType);
         cellCollect3{2,count} = string(fileIndex(:,1:11,d));
         cellCollect3{3,count} = [meanOnMax meanOnMaxOri(1) meanOnNull meanOnNullOri meanOffMax meanOffMaxOri(1) meanOffNull meanOffNullOri];

    
    
        end


         
        holdName = [string(holdName); string(cellType)];
        
        
        
end
%% now use the collected statistics for each parameter set
% cellNames = [cellCollect{1,:}]';
% cellFiles = [cellCollect{2,:}]';
% cellStats = [cellCollect{3,:}]';

% cellNames = [cellCollect2{1,:}]';
% cellFiles = [cellCollect2{2,:}]';
% cellStats = [cellCollect2{3,:}]';
% 
% cellNames = [cellCollect3{1,:}]';
% cellFiles = [cellCollect3{2,:}]';
% cellStats = [cellCollect3{3,:}]';

cellNames = [cellCollect4{1,:}]';
cellFiles = [cellCollect4{2,:}]';
cellStats = [cellCollect4{3,:}]';
% typestoFind = ["On Parasol";"Off Parasol";"On Smooth";"Off Smooth";"B
typesToFind = unique(cellNames);

for z = 1:length(typesToFind)
cellStatsReshaped = reshape(cellStats,8,size(cellStats,1)/8); %this needs to be backwards so you can rotate with ' ... maybe flip(cellStats) would have worked idk
cellStatsReshaped=cellStatsReshaped';

typeName = typesToFind(z);
whichType = strcmp(cellNames,typeName);

statsCollect = cellStatsReshaped(whichType,:);
cellFiles(whichType)
statsForIndex = statsCollect(:,[1 3 5 7]);

filterSmallValues = statsForIndex(:,1)+statsForIndex(:,2);
statsForIndex(filterSmallValues <= 5,1)=0;
statsForIndex(filterSmallValues <= 5,2)=0;


onIndex1 = (statsForIndex(:,1)-statsForIndex(:,2))./(statsForIndex(:,1)+statsForIndex(:,2));
onIndex1(isnan(onIndex1))=0

filterSmallValues = statsForIndex(:,3)+statsForIndex(:,4);
statsForIndex(filterSmallValues <= 5,3)=0;
statsForIndex(filterSmallValues <= 5,4)=0;


offIndex1=(statsForIndex(:,3)-statsForIndex(:,4))./(statsForIndex(:,3)+statsForIndex(:,4));

offIndex1(isnan(offIndex1))=0

onIndexFinal = mean(onIndex1)
onIndexSEM = sem(onIndex1);
offIndexFinal = mean(offIndex1)
offIndexSEM = sem(offIndex1);


% finalCollection(z,1) = onIndexFinal;
% finalCollection(z,2) = offIndexFinal;
% finalError(z,1) = onIndexSEM;
% finalError(z,2) = offIndexSEM;
% finalTypes = typesToFind;


% finalCollection2(z,1) = onIndexFinal;
% finalCollection2(z,2) = offIndexFinal;
% finalError2(z,1) = onIndexSEM;
% finalError2(z,2) = offIndexSEM;
% finalTypes2 = typesToFind;
% 
% finalCollection3(z,1) = onIndexFinal;
% finalCollection3(z,2) = offIndexFinal;
% finalTypes3 = typesToFind;
% finalError3(z,1) = onIndexSEM;
% finalError3(z,2) = offIndexSEM;
% 
% finalCollection4(z,1) = onIndexFinal;
% finalCollection4(z,2) = offIndexFinal;
% finalTypes4 = typesToFind;
finalError4(z,1) = onIndexSEM;
finalError4(z,2) = offIndexSEM;
end

% save('C:\Users\reals\Documents\PhD 2021\ClarinetExports\Oriented Bars\oBarsFinal.mat','finalCollection','finalCollection2','finalCollection3','finalCollection4')
%% now... build the graph?

for w = 1:length(holdName)

onGather = [finalCollection(strcmp(finalTypes,holdName(w)),1) finalCollection2(strcmp(finalTypes2,holdName(w)),1) finalCollection3(strcmp(finalTypes3,holdName(w)),1) finalCollection4(strcmp(finalTypes4,holdName(w)),1)]
offGather = [finalCollection(strcmp(finalTypes,holdName(w)),2) finalCollection2(strcmp(finalTypes2,holdName(w)),2) finalCollection3(strcmp(finalTypes3,holdName(w)),2) finalCollection4(strcmp(finalTypes4,holdName(w)),2)]

if isempty(onGather)
    onGather = 0;
end

if isempty(offGather)
    offGather = 0;
end

finalMean(w,1) = max(onGather);
finalMean(w,2) = max(offGather);

end

%% jk build it now

[one two] = sort(finalMean(:,1));


plot(flip(one));hold on
plot(flip(finalMean(two,2)))
set(gca,'xtick',[1:18],'xticklabel',flip(holdName(two)))




