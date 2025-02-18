%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = getDirectory();
protocolToAnalyze = 'Motion And Noise_Mike';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');
cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
% folderSet(:,:,1:3)=[]; %removes macOS hidden stuff
folderSet(:,:,1:2)=[]; 
%%
 timings = [250,10000,250]; %same for these
 nonlinearityBins = 100;
kmeansGroups=[];
cellnameList = [];
errorGroups = [];

binRate = 1e3;
tFilterMeanData=[];
count=0;
    
    cellTypeData = [];
    typeInd = 1; 
% size(folderSet,3)
% includeTypes = [5 3 2 4 1];
includeTypes = [4];
maxFilterPerCell=[];
minFilterPerCell=[];
timeBetweenPeaksPerCell = [];
peakToPeakMagPerCell = [];
inflectionPointPerCell = [];
startEndDeltaPerCell = [];
initialSlopePerCell = [];
for c = includeTypes
    
    cd(strcat(typesFolder,folderSet(:,:,c))) 
    
    cellType = strtrim(folderSet(:,:,c)); 

    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    
    fileIndex = fileSet(:,:,fileSet(:,1,:) == cellType(1));
    
    hzRSE = [];
    gainRSE = [];
    pcChangeHzSequential= [];
    pcChangeHzRandom = [];
    pcChangeGainSequential = [];
    pcChangeGainRandom = [];
    
    multiCellSeq = [];
    multiCellRand = [];
    multiCellStatic = [];
    
    mseErrorHz = [];
    mseErrorGain = [];
    
    maxFilter = [];
    minFilter = [];
    timeBetweenPeaks= [];
    peakToPeakMag = []; 
    inflectionPoint = [];
    startEndDelta = [];
    initialSlope = [];
    

    for d = 1:size(fileIndex,3)
%         for d = 6:6
           

         load(strtrim(fileIndex(:,:,d)))
         
         stimulusParams.cellName
         %I used to generate filters and NL until I realized Mike had
         %already done them
%   
%             
%             
%             
%             prePts = 250 * 1e-3 * binRate;
%             stimPts = 10000 * 1e-3 * binRate;
%             tailPts = 250 * 1e-3 * binRate;
%            lfilter=[];
% 
%                  stimIndex = prePts + (1 : stimPts);
%        
%                 
%             
%                 R = stimulusParams.response;
% %                
%                 P = zeros(size(R));
%                 for p = 1 : size(R,1)
%                         P(p,:)=stimulusParams.prediction(1,:);
% %                         P(p,:)=P(p,:)/std(P(p,:));
%                 end
%                 
%                 logicalsMike=zeros(3,size(R,1));
%                 
%                 logicalsMike(1,:)=strcmp(stimulusParams.backgroundClasses,'random');
%                 logicalsMike(2,:)=strcmp(stimulusParams.backgroundClasses,'sequential');
%                 logicalsMike(3,:)=strcmp(stimulusParams.backgroundClasses,'stationary');
%                 
%                 if sum(strcmp(stimulusParams.backgroundClasses,'random')) ==0
%                 logicalsMike(1,:)=strcmp(stimulusParams.backgroundClasses,'jittering');
%                 logicalsMike(2,:)=strcmp(stimulusParams.backgroundClasses,'drifting');
%                 logicalsMike(3,:)=strcmp(stimulusParams.backgroundClasses,'stationary');
%                 end
%                 
%                 
%                for g = 1:3 
%                 [xBin(g,:),yBin(g,:)] = binNonlinearity(P(logical(logicalsMike(g,:)),stimIndex(501:end)),R(logical(logicalsMike(g,:)),stimIndex(501:end)),nonlinearityBins);
%                end
%                 
% %                 nlParams(g,:) = models.ln.fitNonlinearityParams(xBin(g,:), yBin(g,:));
%                 prediction(indices,:) = P;
% 

% silly if/else statement but it's what I used
if strcmp(analysisParams.classLabels(1),'random') %most files have random as row 1 because of how stimulus is written.  I kept that format in analysis     
    xBin(1,:) = analysisParams.xBin(1,:);
    yBin(1,:) = analysisParams.yBin(1,:);
    xBin(2,:)=analysisParams.xBin(2,:);
    yBin(2,:)=analysisParams.yBin(2,:);
    xBin(3,:)=analysisParams.xBin(3,:);
    yBin(3,:)=analysisParams.yBin(3,:);
else %this puts the random condition as first row in array for the few cases where motion is first row (b/c of early stim versions) Ultimately, 1 = random 2 = motion 3 = static
    xBin(2,:) = analysisParams.xBin(1,:);
    yBin(2,:) = analysisParams.yBin(1,:);
    xBin(1,:)=analysisParams.xBin(2,:);
    yBin(1,:)=analysisParams.yBin(2,:);
    xBin(3,:)=analysisParams.xBin(3,:);
    yBin(3,:)=analysisParams.yBin(3,:); 
end



% % % %% for TF gazing

figure(12)
subplot(2,ceil(size(fileIndex,3)/2),d)
plot(analysisParams.lfilter(1:300))

figure(14)
hold on
plot(analysisParams.lfilter(1:100))

lfilter=analysisParams.lfilter(1:500);

indexForEnd = find(lfilter(:).*circshift(lfilter(:), [-1 0]) <= 0);
indexForEnd = indexForEnd(indexForEnd >200);
indexForStart = find(lfilter(:).*circshift(lfilter(:), [-1 0]) <= 0);
indexForStart = indexForStart(indexForStart<20);

if isempty(indexForStart)
    indexForStart = 1;
end

if isempty(indexForEnd)
    indexForEnd = 0;
end


maxFilter(d) = max(lfilter);
minFilter(d) = min(lfilter);
timeBetweenPeaks(d) = abs(find(lfilter==maxFilter(d)) - find(lfilter==minFilter(d)));
peakToPeakMag(d) = maxFilter(d) + abs(minFilter(d)); 
startEndDelta(d) = indexForEnd(1) - indexForStart(end);
if startEndDelta(d) < 0
    startEndDelta(d) = [];
end
initialSlope(d) = ((lfilter(find(diff(lfilter)==max(diff(lfilter)))+1)-lfilter(find(diff(lfilter)==max(diff(lfilter)))-1)))/2; %rise/run


[S, L] = bounds(diff(lfilter));
if find(diff(lfilter)==S) > find(diff(lfilter)==L)
    inflectionPoint(d)= lfilter(find(diff(lfilter)==S)); 
else
    inflectionPoint(d)= lfilter(find(diff(lfilter)==L));
end


% % %% for NL debuggin'

  figure(13)
  colororder({'b','r','k'})
  subplot(2,ceil(size(fileIndex,3)/2),d)
  plot(xBin',yBin')
  
% 

if d == 16
    figure(33)
    plot(xBin',yBin')
    figure(34)
    plot(lfilter)
    
    for e = 1:(size(xBin,1))
    nlParamsExample(e,:) = fitNonlinearityParams(xBin(e,:), yBin(e,:));
 
    figure(35);hold on;plot(xBin(e,:),outputNonlinearity(nlParamsExample(e,:),xBin(e,:)),'--') 
    
    
    end
    
end
%%%%%% get spike rates

    if iscell(stimulusParams.response)
        meanRand = mean(sum(stimulusParams.response{2,1},2));
        meanSeq = mean(sum(stimulusParams.response{2,2},2));
        meanStatic = mean(sum(stimulusParams.response{2,3},2));
    else
        allSpikes = stimulusParams.response(:,250:end);
        bgClasses = stimulusParams.backgroundClasses;
        seqSpikes = allSpikes(contains(bgClasses,["sequential","drifting"]),:);
        
        randSpikes = allSpikes(contains(bgClasses,["random","jittering"]),:);
        staticSpikes = allSpikes(contains(bgClasses,"stationary"),:);

        meanSeq = mean(sum(seqSpikes,2));
        meanRand = mean(sum(randSpikes,2));
        meanStatic = mean(sum(staticSpikes,2));
    end
multiCellSeq(d) = meanSeq;
multiCellRand(d) = meanRand;
multiCellStatic(d) = meanStatic;


%%%%%
                %go with 1 being random, 2 sequential, 3 static
              starterParams = [.2,-1.5];
           

                %All
               
              %each half modeled for broad thorny, only Off half used
              if strcmp(cellType,'Broad Thorny')
                 mParams = fitMultiVarParams(xBin([1,2,3],51:100),yBin([1,2,3],51:100),1,starterParams);
                 outNLHZOn = multiHZNL(mParams,xBin([1,2,3],51:100));

%                  figure;colororder(['b';'r';'k']);plot(xBin(:,51:100)',outNLHZOn','--');hold on;plot(xBin(:,51:100)',yBin(:,51:100)')
                 

                 mParamsHz = fitMultiVarParams(xBin([1,2,3],1:50),yBin([1,2,3],1:50),1,starterParams);
                 outNLHZOff = multiHZNL(mParamsHz,xBin([1,2,3],1:50));
                 
                 mParamsGain = fitMultiVarParams(xBin([1,2,3],1:50),yBin([1,2,3],1:50),0,starterParams);
                 outNLGainOff = multiGainNL(mParamsGain,xBin([1,2,3],1:50));
                 
                 

%                  colororder(['b';'r';'k']);plot(xBin(:,1:50)',outNLHZOff','--');hold on;plot(xBin(:,1:50)',yBin(:,1:50)')
%                  title('Broad Thorny NLs w/ Separate On/Off Model (XshiftOnly)')
%                  xlabel('Spike Rate (Hz)')
%                  ylabel('Input')
%                  
                 
%                  figure;colororder(['b';'r';'k']);plot(xBin(:,1:50)',outNLGainOff','--');hold on;plot(xBin(:,1:50)',yBin(:,1:50)')
%                  title('Broad Thorny NLs w/ Separate On/Off Model (XshiftOnly)')
%                  xlabel('Spike Rate (Hz)')
%                  ylabel('Input')
                
                mseErrorHz(d) = immse(yBin(2,1:50),outNLHZOff(2,:));
                mseErrorGain(d) = immse(yBin(2,1:50),outNLGainOff(2,:));
                
                figure(97);clf
                plot(xBin(2,1:50),yBin(2,1:50))
                hold on
                plot(xBin(2,1:50),outNLHZOff(2,:),'--')
                figure(96);clf
                plot(xBin(2,1:50),yBin(2,1:50))
                hold on
                plot(xBin(2,1:50),outNLGainOff(2,:))
                
                 
              else

              mParamsHz = fitMultiVarParams(xBin([1,2,3],:),yBin([1,2,3],:),1,starterParams);
              outNLHZ = multiHZNL(mParamsHz,xBin([1,2,3],:));
              
              
              
              mParamsGain = fitMultiVarParams(xBin([1,2,3],:),yBin([1,2,3],:),0,starterParams);
              outNLGain = multiGainNL(mParamsGain,xBin([1,2,3],:));
              
              if d==6
                  figure(74);colororder(['b';'r';'k']);plot(xBin',outNLHZ','--',xBin',yBin');
                  figure(75);colororder(['b';'r';'k']);plot(xBin',outNLGain','--',xBin',yBin');
              end    
              
              mseErrorHz(d) = immse(yBin(2,:),outNLHZ(2,:));
              mseErrorGain(d) = immse(yBin(2,:),outNLGain(2,:));
              
%                figure(97);clf 
%                 plot(xBin(2,1:50),yBin(2,1:50))
%                 hold on
%                 plot(xBin(2,1:50),outNLHZOff(2,:),'--')
%                 figure(96);clf
%                 plot(xBin(2,1:50),yBin(2,1:50))
%                 hold on
%                 plot(xBin(2,1:50),outNLGainOff(2,:))
%                 pause
                
              end
              
              hzRSE(d,:) = [mParamsHz(3),mParamsHz(4),mParamsHz(5)];
              gainRSE(d,:) = [mParamsGain(3),mParamsGain(4),mParamsGain(5)];
              
              
              

              
%               line([0 1],[0 1],'Color','k')
%               
              
%               figure(19)
%               plot(xBin(1,:),outNLHZ(1,:),'b')
%               hold on
%               plot(xBin(2,:),outNLHZ(2,:),'r')
%          
%               plot(xBin(1,:),yBin(1,:),'b--')
%               plot(xBin(2,:),yBin(2,:),'r--')
%               legend('Random-Model','Sequential-Model','Random-Data','Sequential-Data')
%               title('Raw and Modeled NLs - XShift')
%               xlabel('input')
%               ylabel('spike rate (Hz)')




                
        end
%    
%        figure(29)
%        pcChangeHzShift =  (-diff(hzRSE,1,2)./sum(hzRSE,2))*100;
%        plot(ones(size(pcChangeHzShift)),pcChangeHzShift,'.','LineWidth',10)
%        xlabel('arbitrary')
%        ylabel('Percent Change Horizontal Shift Parameter')
%        title('Change in X Shift Parameter in NL Models for Motion & Random Conditions')

%%% temporal filter means

%maxFilter,minFilter,timeBetweenPeaks,peakToPeakMag



maxFilterPerCell = [maxFilterPerCell mean(maxFilter)];
minFilterPerCell = [minFilterPerCell mean(minFilter)];
timeBetweenPeaksPerCell = [timeBetweenPeaksPerCell mean(timeBetweenPeaks)];
peakToPeakMagPerCell = [peakToPeakMagPerCell mean(peakToPeakMag)];
inflectionPointPerCell = [inflectionPointPerCell mean(inflectionPoint)];
startEndDeltaPerCell = [startEndDeltaPerCell mean(startEndDelta)];
initialSlopePerCell = [initialSlopePerCell mean(initialSlope)];





count = count+1;

seqChange=[];
randChange=[];
seqChange = 100*(multiCellSeq-multiCellStatic)./multiCellStatic;
randChange = 100*(multiCellRand-multiCellStatic)./multiCellStatic;

seqX = ones(1,length(multiCellSeq))*count;
figure(11);hold on
% seqX=seqX+1;
changeY = [seqChange' randChange'];
changeX = [seqX' seqX'+.2];
plot(changeX, changeY,'.')
axis([.5 5.5 -50 120])
line([count' count'+.2],[seqChange(1:end)' randChange(1:end)'])
line([count' count'+.2],[mean(seqChange) mean(randChange)],'LineWidth',3)
% makeAxisStruct(gca,strtrim(['changeDataMikeLINES' 'offSmooth']))

errorSequential(c) = sem(seqChange);
errorRandom(c) = sem(randChange);
%         cellType = fileIndex(:,1:10,1);



%%%%%%%%%%%%%%%%%%% UNITY LINE GRAPHS
mseErrorHzMean = mean(mseErrorHz);
mseErrorGainMean = mean(mseErrorGain);

mseErrorHzError(d) = sem(mseErrorHz);
mseErrorGainError(d) = sem(mseErrorGain);

colorSet = ['b.'; 'r.'];
colorSet2 = ['b*'; 'r*'];

figure(98);hold on
dotSize = 500;
scatter(mseErrorHz,mseErrorGain,dotSize,'k.')
% line([0 20],[0 20],'Color','k')
axis([0 25 0 25])

figure(92);hold on
scatter(mseErrorHzMean,mseErrorGainMean,dotSize*3,'k+')
line([0 10],[0 10],'Color','k')
axis([0 10 0 10])
% meanErrorHZ = mean(MSEseqHZ);
% semErrorHZ=sem(MSEseqHZ);
% 
% meanErrorGain = mean(MSEseqGain);
% semErrorGain = sem(MSEseqGain);
% 
% 
% figure
% changeY = [1-(MSEseqHZ/max([MSEseqHZ MSEseqGain]))'];
% changeX = [1-(MSEseqGain/max([MSEseqHZ MSEseqGain]))'];
% plot(changeX, changeY,'.')
% line([0 1],[0 1],'Color','k')
% title(cellType)
% hold on
% changeY = [1-MSErandHZ/max([MSErandHZ MSErandGain])' ];
% changeX = [1-MSErandGain/max([MSErandHZ MSErandGain])'];
% plot(changeX, changeY,'.')
% 
% figure
% changeY = [MSEseqHZ'];
% changeX = [MSEseqGain'];
% plot(changeX, changeY,'.')
% hold on
% changeY = [MSErandHZ'];
% changeX = [MSErandGain'];
% plot(changeX, changeY,'.')
% line([0 max([MSErandGain MSErandHZ MSEseqHZ MSEseqGain])],[0 max([MSErandGain MSErandHZ MSEseqHZ MSEseqGain])],'Color','k')
% title(strtrim([cellType ' Straight MSE']))
% plot(meanErrorGain,meanErrorHZ)
% line([(meanErrorGain-semErrorGain) (meanErrorGain+semErrorGain)], [meanErrorHZ meanErrorHZ])
% line([(meanErrorGain) (meanErrorGain)], [meanErrorHZ-semErrorHZ meanErrorHZ+semErrorHZ])
% % makeAxisStruct(gca,strtrim(['MikeUnity_' 'onSmoothNorm']))
% 
% meanMSESeqHZ(typeInd)=mean(MSEseqHZ);
% semMSESeqHZ(typeInd)=sem(MSEseqHZ);
% meanMSESeqGain(typeInd)=mean(MSEseqGain);
% semMSESeqGain(typeInd)=sem(MSEseqGain);
% meanMSErandHZ(typeInd)=mean(MSErandHZ);
% meanMSErandGain(typeInd)=mean(MSErandGain);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%graphing
      figure(29) % seq static
       pcChangeHzSequential =  (diff(hzRSE(:,[2 3]),1,2)./sum(hzRSE(:,[2 3]),2))*100;
       plot(count,mean(pcChangeHzSequential),'r.','MarkerSize',30)
       hold on
       plot(count,pcChangeHzSequential,'ro', 'MarkerSize',5)
       

      errorHzSeq(c) = sem(pcChangeHzSequential);
      
      
    
        


       pcChangeHzRandom =  (diff(hzRSE(:,[1 3]),1,2)./sum(hzRSE(:,[1 3]),2))*100;
       plot(count+.2,mean(pcChangeHzRandom),'b.','MarkerSize',30)
       plot(count+.2,pcChangeHzRandom,'bo', 'MarkerSize',5)
       xlabel('cell type')
       ylabel('Percent Change')
       title('Extent of XShift Parameter Change Relative to Static Condition')
       

    
errorHzRand(c) = sem(pcChangeHzRandom);
    
       
      figure(31) %seq static
       pcChangeGainSequential =  (diff(gainRSE(:,[3 2]),1,2)./sum(gainRSE(:,[3 2]),2))*100;
       plot(count,mean(pcChangeGainSequential),'r.','MarkerSize',30)
       hold on
       plot(count,pcChangeGainSequential,'ro','MarkerSize',5)
       
       
errorGainSeq(c) = sem(pcChangeGainSequential);
       

       pcChangeGainRandom =  (diff(gainRSE(:,[3 1]),1,2)./sum(gainRSE(:,[3 1]),2))*100;
       plot(count+.2,mean(pcChangeGainRandom),'b.','MarkerSize',30)
       plot(count+.2,pcChangeGainRandom,'bo', 'MarkerSize',5)
       xlabel('cell type')
       ylabel('Percent Change')
       title('Extent of Gain Parameter Change Relative to Static Condition')

           
  errorGainRand(c) = sem(pcChangeGainRandom);

end

tFilterMeanData{1,1} = 'Filter Metric';
tFilterMeanData{2,1} = 'Metric Mean Per Type';
tFilterMeanData{1,2} = 'Max Filter';
tFilterMeanData{2,2} = maxFilterPerCell;
tFilterMeanData{1,3} = 'Min Filter';
tFilterMeanData{2,3} = minFilterPerCell;
tFilterMeanData{1,4} = 'Time Between Peaks';
tFilterMeanData{2,4} = timeBetweenPeaksPerCell;
tFilterMeanData{1,5} = 'Peak to Peak Mag';
tFilterMeanData{2,5} = peakToPeakMagPerCell;
tFilterMeanData{1,6} = 'Start to End Delta';
tFilterMeanData{2,6} = startEndDeltaPerCell;
tFilterMeanData{1,7} = 'Initial Slope';
tFilterMeanData{2,7} = initialSlopePerCell; %broken



% 
% tFilterMeanData{1,c+1} 
% 
% tFilterMeanData{1,c+1} = cellType;
% tFilterMeanData{2,c+1} = [mean(maxFilter) mean(minFilter) mean(timeBetweenPeaks) mean(peakToPeakMag) mean(inflectionPoint)];
% tFilterMeanData{3,c+1} = [sem(maxFilter) sem(minFilter) sem(timeBetweenPeaks) sem(peakToPeakMag) sem(inflectionPoint)];
% 
% tFilterMeanData{4,2} = ['Filter Max ', 'Filter Min ',];
% tFilterMeanData{4,3} = [' Time Between Peaks', ' Peak to Peak Mag'];
% tFilterMeanData{4,4} = ['Inflection Point'];
% 
% maxFilterMeanCollection = [tFilterMeanData{2,2:end}];
% 
% minFilterMeanCollection = maxFilterMeanCollection(1,[2 7 12 17 22]);
% deltaPeaks = maxFilterMeanCollection(1,[3 8 13 18 23]);
% peaktoPeakMag = maxFilterMeanCollection(1,[4 9 14 19 24]);
% inflectionPointMean = maxFilterMeanCollection(1,[5 10 15 20 25]);
% maxFilterMeanCollection = maxFilterMeanCollection(1,[1 6 11 16 21]);



figure(29)
hzFig = gca;
hzFig.XTick =[1.2 2.2 3.2 4.2 5.2];
xticklabels({'OnSmooth','OnParasol','OffParasol','OffSmooth','BT'})
figure(31)
gainFig = gca;
gainFig.XTick =[1.2 2.2 3.2 4.2 5.2];
xticklabels({'OnSmooth','OnParasol','OffParasol','OffSmooth','BT'})
% X=[];
% Y=[];
% for g= 1:3
% nlParams(g,:) = models.ln.fitNonlinearityParams(X(g,:), Y(g,:));
% plot(X(g,:),outputNonlinearity(nlParams(g,:),X(g,:)))
% hold on
% end
% X=X'
% Y=Y'
%%
save('MeansandErrorsofError_Mike','meanMSESeqHZ','semMSESeqHZ','meanMSESeqGain','semMSESeqGain')
%% Scatter!
colors = [1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;.5 .2 .3;.7 .7 .7;.1 .2 .3;.4 .1 .8;.2 .3 .6;0 .1 .9;.4 .9 .1];
seqX = ones(1,length(multiCellSeq));
randX = 2*ones(1,length(multiCellRand));
staticX = 3*ones(1,length(multiCellStatic));

scatter([seqX randX staticX],[multiCellSeq multiCellRand multiCellStatic])
hold on
scatter([1 2 3],[mean(multiCellSeq) mean(multiCellRand) mean(multiCellStatic)],'SizeData',20,'LineWidth',4)

figure
for j = 1:length(multiCellSeq)
    scatter([1 2 3],[multiCellSeq(j) multiCellRand(j) multiCellStatic(j)],'MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:))
    hold on
end
scatter([1 2 3],[mean(multiCellSeq) mean(multiCellRand) mean(multiCellStatic)],'d','SizeData',20,'LineWidth',4)
%% testing ground

for j = 1:size(staticSpikes,1)
    plot(seqSpikes(j,:))
    pause
end

%% graphing (old)
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

%% mistake eraser

%     cd(strcat(exportFolder,fileIndex(:,1:4,d),'_',fileIndex(:,5:8,d)))
    splitFactors = ["noiseClass","epochGroupLabel","frameDwell","backgroundClass"];
    load(strcat(fileIndex(:,1:11,d),'.mat'))
    spikeEpochs=epochs;
    load(strcat(fileIndex(:,1:11,d),'_FT.mat'))
    frameTs = epochs;
    saveLabel = 'MN';
    [splitCell,indexHolder,spikingData,frameTimings,metaData] = makeData(spikeEpochs,frameTs,protocolToAnalyze,splitFactors,6,recordingType);
    
    cd(typeDir(1).folder)
    save(strcat(fileIndex(:,1:11,d),'_','MN','.mat'),'splitCell','indexHolder','frameTimings','-append')
