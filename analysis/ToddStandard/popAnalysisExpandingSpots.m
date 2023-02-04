    %popAnalysisExpandingSpots
%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'Expanding Spots';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

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
holdName=[];
expSpotsOut=cell(0);
count=0; 
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
collectedMeans = [];
        for d = 1:size(fileIndex,3)
      count=count+1;

         load(strtrim(fileIndex(:,:,d)))
         
     
         
         % get timings for each cell just in case diff numbers used 
         timings = [metaData.preTime,metaData.stimTime,metaData.tailTime];
         stimStart = (timings(1)*1e-3)*sampleRate+1;
         stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
         stimOff = (timings(1)+timings(2)+timings(3))*10;
         % timings done
         
         sumsort=[];
         clear labelHolder
             for s = 1:size(indexHolder,2)
                 labelHolder(s,:) = indexHolder{1,s};
             end

        
             
 
         spotIntensityOnly = double(labelHolder(:,1));
         bgIntensityOnly = double(labelHolder(:,2));
         posContrast = spotIntensityOnly>=.5;
         negContrast = spotIntensityOnly<.5;
         
         bgOn = bgIntensityOnly==.5;
         bgOff = bgIntensityOnly==0;
        

         
         indPoswithBG= vertcat(indexHolder{2,bgOn&posContrast});
         indNegwithBG = vertcat(indexHolder{2,bgOn&negContrast});
         indNoBG = vertcat(indexHolder{2,bgOff&posContrast});
         
         indPoswithBGRadii = double(splitCell{2,3}(indPoswithBG));
         indNegwithBGRadii = double(splitCell{2,3}(indNegwithBG));
         indNoBGRadii = double(splitCell{2,3}(indNoBG));
         
         radii=unique(double(splitCell{2,3}));
         radiiMean = zeros(size(radii));
         for xx = 1:length(radii)
            
            if length(indPoswithBG)>length(indNegwithBG)
            arrayIndex= indPoswithBG(indPoswithBGRadii==(radii(xx)),1);
            radiiMean(xx) = mean(sum(spikingData(arrayIndex,stimStart:stimOff),2));
            elseif length(indNoBG) > length(indPoswithBG) && length(indNoBG) > length(indNegwithBG)
            arrayIndex= indNoBG(indNoBGRadii==(radii(xx)),1);
            radiiMean(xx) = mean(sum(spikingData(arrayIndex,stimStart:stimOff),2));
            else
            arrayIndex= indNegwithBG(indNegwithBGRadii==(radii(xx)),1);
            radiiMean(xx) = mean(sum(spikingData(arrayIndex,stimStart:stimOff),2));
            end

         end
         
         if sum(radiiMean)==0 
             count = count -1;
             continue   
         end
         
         if sum(isnan(radiiMean))>0
             count=count-1;
             continue
         end
             
         
         afterMax = radiiMean(find(radiiMean==max(radiiMean)):end);
         afterMaxRadii = radii(find(radiiMean==max(radiiMean)):end);
         suppressTest = mean(radiiMean((length(radiiMean)-2):end));
         suppressionIndex = (max(radiiMean)-suppressTest)/max(radiiMean);
         percentReduction = 1 - min(afterMax)/max(radiiMean);
         maxRadius = min(radii(radiiMean==max(radiiMean)));
         
         halfMaxRadius = radii(round(find(radiiMean==max(radiiMean))/2));
         riseRateLate = (radiiMean(find(radiiMean==max(radiiMean))) - radiiMean(round(find(radiiMean==max(radiiMean))/2)))/ (maxRadius-halfMaxRadius);
         riseRateEarly = (radiiMean(round(find(radiiMean==max(radiiMean))/2)) - radiiMean(1)) / (halfMaxRadius-radii(1));
         suppressRate = (radiiMean(find(radiiMean==max(radiiMean))) - min(afterMax))/ (maxRadius - min(afterMaxRadii(afterMax == min(afterMax))));
         
         
         expSpotsOut{1,count} = string(cellType);
         expSpotsOut{2,count} = string(fileIndex(:,1:11,d));
         expSpotsOut{3,count} = [maxRadius(1), halfMaxRadius(1), percentReduction(1), riseRateEarly(1), riseRateLate(1), suppressRate(1) suppressionIndex(1)];
         
         collectedMeans(d,1:length(radiiMean)) = radiiMean;

        holdFile = [string(holdFile) string(fileIndex(:,1:11,d))];
        end
        
        meanType = mean(collectedMeans,1);
        areaCurve(c,1:length(meanType))=meanType;
        
        %must clear zeros from above because matlab puts zeros to
        %compensate for changing size of array #thanksmatlab (actually)
        areaCurve(c,length(meanType):size(areaCurve,2))=areaCurve(c,length(meanType)); 
        
        areaRadii(c,1:length(radii))=radii;
        areaRadii(c,areaRadii(c,:)==0)=max(areaRadii(c,:));
       
         
         
        holdName = [string(holdName); string(cellType)];
        
        
        
        
end
% 
 cd(exportFolder)

%   plot(areaRadii',areaCurve(:,1:size(areaRadii,2))'./max(areaCurve(:,1:size(areaRadii,2))'));legend(holdName)
plot(areaRadii(1:8,:)',areaCurve(1:8,1:size(areaRadii,2))'./max(areaCurve(1:8,1:size(areaRadii,2))'));legend(holdName(1:8))
% cellF1F2 = normr(cellF1F2); 
% cellF1 = normr(cellF1); 
% cellF2 = normr(cellF2); 


% save(strcat('smtfF1F2.mat'),'cellF1F2','holdName')
% save(strcat('smtfF1.mat'),'cellF1','holdName')
% save(strcat('smtfF2.mat'),'cellF2','holdName')

save('Expanding Spots Stats','expSpotsOut')
%% plotter
%choose cell types to plot
cellTypeLogical = contains(holdName,["ON Parasol","Off Parasol","On Smooth","Off Smooth"]);
plot(areaRadii(cellTypeLogical,:)',areaCurve(cellTypeLogical,1:size(areaRadii,2))');legend(holdName(cellTypeLogical))
%% check stats

typeOut = vertcat(expSpotsOut{1,:});
fileOut = vertcat(expSpotsOut{2,:});
statsOut = vertcat(expSpotsOut{3,:});
% "Off Parasol","On Smooth","Off Smooth"

typeOutLogical = contains(typeOut,["strong on OS"]);
statsOut(isnan(statsOut(:,1:6)))=0;

tsneExpSpots = tsne(statsOut(:,3:4));

gscatter(tsneExpSpots(~strcmp(typeOut,"Unknown"),1),tsneExpSpots(~strcmp(typeOut,"Unknown"),2),typeOut(~strcmp(typeOut,"Unknown")),colors2)

% typeLogical = contains(typeOut,"ON Parasol");

% plot(statsOut(typeLogical,3:5)')

%% Load smtf spot and expanding spot

cd('/Users/reals/Documents/PhD 2021/ClarinetExports/S MT Fspot/')
load('sMTF Spot Size.mat')
cd('/Users/reals/Documents/PhD 2021/ClarinetExports/Expanding Spots/')
load('Expanding Spots Stats.mat')

%% condense cell

allCellTypes = vertcat([sMTFSpotsOut{1,:}]',[expSpotsOut{1,:}]');
allFileNames = vertcat([sMTFSpotsOut{2,:}]',[expSpotsOut{2,:}]');

for cellNumStatsExpSpots = 1:length(expSpotsOut)
    comboStatsES(cellNumStatsExpSpots,:)=expSpotsOut{3,cellNumStatsExpSpots};
end
comboStatsES(isnan(comboStatsES))=0;
count = 0;
for cellNumStatsSMTF = 1:length(sMTFSpotsOut)
    
    if ~isempty(sMTFSpotsOut{3,cellNumStatsSMTF})
        
    count=count+1;    
    comboStatsSMTF(count,:)=sMTFSpotsOut{3,cellNumStatsSMTF};
    end
end
% comboStatsSMTF(isnan(comboStatsSMTF))=0;

allStats = vertcat(comboStatsSMTF,comboStatsES);

cd('/Users/reals/Documents/PhD 2021/ClarinetExports/')
save('combinedProtocolStats','allCellTypes','allFileNames','allStats')
%% combine all

cd('/Users/reals/Documents/PhD 2021/ClarinetExports/')
load('PercentChangeMotionandNoise.mat')
allCellTypesMN=allCellTypes;
allFileNamesMN=string(allCellNames2);
cd('/Users/reals/Documents/PhD 2021/ClarinetExports/')
load('combinedProtocolStats')

allCellTypes(16)=[];
allFileNames(16)=[];
allStats(16,:)=[];
%% match

matchLogical = contains(allFileNames,allFileNamesMN);
matchedFileNames = allFileNames(matchLogical);
% fileSelect=allFileNamesMN(contains(allFileNamesMN,matchedFileNames));

matchedStats = allStats(matchLogical,3);
matchedTypes = allCellTypes(matchLogical);

suppIndex = allIndex(contains(allFileNamesMN,matchedFileNames));

changeDiff=allSeqMeans-allRandMeans;
changeSelect=changeDiff(contains(allFileNamesMN,matchedFileNames));
justSeq = allSeqMeans(contains(allFileNamesMN,matchedFileNames));

typeSelect=allCellTypesMN(contains(allFileNamesMN,matchedFileNames));
fileSelect=allFileNamesMN(contains(allFileNamesMN,matchedFileNames));



[sortedFSMN MNFileIndex] = sort(fileSelect);
[sortedFSAll AllFileIndex] = sort(matchedFileNames);

diffChange = changeSelect(MNFileIndex);
pSeqChange = justSeq(MNFileIndex);
pSuppress = matchedStats(AllFileIndex);
typeFinal = typeSelect(MNFileIndex);
suppressionIndex = matchedStats(AllFileIndex);
adaptIndex = suppIndex(MNFileIndex);

figure
% plot(diffChange(strcmp(typeFinal,"ONParasol"))/max(diffChange(strcmp(typeFinal,"ONParasol"))))
% hold on
% plot(pSuppress(strcmp(typeFinal,"ONParasol")))

% plot(pSeqChange(strcmp(typeFinal,"ONParasol"))/max(pSeqChange(strcmp(typeFinal,"ONParasol"))))
% hold on
% plot(pSuppress(strcmp(typeFinal,"ONParasol")))

% seqX = ones(1,length(pSeqChange))+2;
% seqX = seqX(strcmp(typeFinal,"OFFSmooth"));
% plotDots = [pSeqChange(strcmp(typeFinal,"OFFSmooth"))/max(abs(pSeqChange(strcmp(typeFinal,"OFFSmooth")))) pSuppress(strcmp(typeFinal,"OFFSmooth"))/max(pSuppress(strcmp(typeFinal,"OFFSmooth")))];
% plotDots = [diffChange(strcmp(typeFinal,"ONParasol"))/max(diffChange(strcmp(typeFinal,"ONParasol"))) pSuppress(strcmp(typeFinal,"ONParasol"))/max(pSuppress(strcmp(typeFinal,"ONParasol")))];

% plot([seqX' seqX'+.2],plotDots,'.')
% hold on
% plot(seqX+.2, pSuppress(strcmp(typeFinal,"ONParasol"))/max(pSuppress(strcmp(typeFinal,"ONParasol"))),'.')
% axis([.5 1.7 -1 2])

% line([seqX(1)' seqX(1)'+.2],plotDots)
% makeAxisStruct(gca,strtrim(['linesSuppressionComparison' '_OFFSmooth']))

load('C:\Users\reals\Documents\PhD 2021\ClarinetExports\Doves Movie\colors.mat')
gscatter(suppressionIndex(~strcmp(typeFinal,"Unknown")),pSeqChange(~strcmp(typeFinal,"Unknown")),typeFinal(~strcmp(typeFinal,"Unknown")),colors2)
%%
whatCells=["ONParasol" "OFFParasol" "ONSmooth" "OFFSmooth" "BroadThorny" "Strong ON OS" "Recursive" "lowTFNOTBT"];
gscatter(suppressionIndex(contains(typeFinal,whatCells)),pSeqChange(contains(typeFinal,whatCells)),typeFinal(contains(typeFinal,whatCells)),colors2)
makeAxisStruct(gca,strtrim(['suppressionVSfacilitation']))
%% 
whichMean = "OFFSmooth";
useLogical = strcmp(typeFinal,whichMean);
mean(suppressionIndex(useLogical))
sem(suppressionIndex(useLogical))
mean(pSeqChange(useLogical))
sem(pSeqChange(useLogical))
%%
tsneExpSpots = tsne(comboStats(:,1:2));

gscatter(tsneExpSpots(~strcmp(types1,"Unknown"),1),tsneExpSpots(~strcmp(types1,"Unknown"),2),types1(~strcmp(types1,"Unknown"))',colors2)
%%

%%
% '[maxRadius(1), halfMaxRadius(1), percentReduction(1), riseRateEarly(1), riseRateLate(1), suppressRate(1)];'
typeWanted = 'Off Smooth';
typeList = [expSpotsOut{1,:}]';
dataList = vertcat(expSpotsOut{3,:});
typeL = strcmp(typeList,typeWanted);
onlyTypeData = dataList(typeL,:);
histogram(onlyTypeData(:,1),'BinWidth',25)
xlabel('Max Response Radius')
ylabel('number')
title(strcat(typeWanted,'- distribution of area curve peaks'))
% figure
% histogram(onlyTypeData(:,2),'BinWidth',25)
figure
scatter(ones(size(onlyTypeData(:,3))),onlyTypeData(:,3))
title(strcat(typeWanted,' - Area Summation Extent of Suppression'))
label = {'Percent Reduction'};
set(gca,'xtick',[1],'xticklabel',label)
ylabel('Percent Reduction')
figure
scatter(ones(size(onlyTypeData(:,3))),onlyTypeData(:,4))
hold on
scatter(ones(size(onlyTypeData(:,3)))+1,onlyTypeData(:,5))
scatter(ones(size(onlyTypeData(:,3)))+2,abs(onlyTypeData(:,6)))
axis([0 4 0 .55])
names = {'Rise Rate - Early'; 'Rise Rate - Late'; 'SuppressRate'};
set(gca,'xtick',[1:3],'xticklabel',names)
title(strcat(typeWanted,' - Area Summation Rise and Suppression'))
ylabel('slope')