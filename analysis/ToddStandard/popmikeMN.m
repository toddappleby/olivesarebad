%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
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

    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == 'P');
    
    cellTypeData = [];
    typeInd = 1;
%     

        for d = 1:size(fileIndex,3)
           

         load(strtrim(fileIndex(:,:,d)))
         
         stimulusParams.cellName
         %Where params are most appropriate?  prob the ones I ran the most:
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
if strcmp(analysisParams.classLabels(1),'random')
            
xBin(1,:) = analysisParams.xBin(1,:);
yBin(1,:) = analysisParams.yBin(1,:);
xBin(2,:)=analysisParams.xBin(2,:);
yBin(2,:)=analysisParams.yBin(2,:);
xBin(3,:)=analysisParams.xBin(3,:);
yBin(3,:)=analysisParams.yBin(3,:);
else
xBin(2,:) = analysisParams.xBin(1,:);
yBin(2,:) = analysisParams.yBin(1,:);
xBin(1,:)=analysisParams.xBin(2,:);
yBin(1,:)=analysisParams.yBin(2,:);
xBin(3,:)=analysisParams.xBin(3,:);
yBin(3,:)=analysisParams.yBin(3,:);
    
end
            
% yBin(1,xBin(1,:)<0)=0;
% yBin(2,xBin(2,:)<0)=0;
% yBin(3,xBin(3,:)<0)=0;

%%%%%% spike count stuff -- actuallty spike rate counts because psth

allSpikes = stimulusParams.response(:,500:end);
bgClasses = stimulusParams.backgroundClasses;
seqSpikes = allSpikes(contains(bgClasses,["sequential","drifting"]),:);
randSpikes = allSpikes(contains(bgClasses,["random","jittering"]),:);
staticSpikes = allSpikes(contains(bgClasses,"stationary"),:);

meanSeq = mean(sum(seqSpikes,2));
meanRand = mean(sum(randSpikes,2));
meanStatic = mean(sum(staticSpikes,2));

multiCellSeq(d) = meanSeq;
multiCellRand(d) = meanRand;
multiCellStatic(d) = meanStatic;

%%%%%
                %go with 1 being random, 2 sequential, 3 static
              starterParams = [.2,-1.5];
              %horizontal  
%                 % rand seq
%               mParams = fitMultiVarParams(xBin([1,2],:),yBin([1,2],:),1,starterParams);
%               outNLHZ = multiHZNL(mParams,xBin([1,2],:));
%               hzRSE = [mParams(3),mParams(4)];
%               
                % rand static
              mParams = fitMultiVarParams(xBin([1,3],:)',yBin([1,3],:)',1,starterParams);
              outNLHZ = multiHZNL(mParams,xBin([1,3],:)');
              staticErrorHZ1(d) = immse(yBin(3,:)',outNLHZ(:,2));
              
           
            MSErandHZ2(d) = immse(yBin(1,:)',outNLHZ(:,1));
            MSErandHZ(d) = MSErandHZ2(d)/staticErrorHZ1(d);
             
              
                % seq static
              mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',1,starterParams);
              outNLHZ = multiHZNL(mParams,xBin([2,3],:)');
              

              staticErrorHZ2(d) = immse(yBin(3,:)',outNLHZ(:,2));
            MSEseqHZ2(d) = immse(yBin(2,:)',outNLHZ(:,1));
            MSEseqHZ(d) = MSEseqHZ2(d)/staticErrorHZ2(d);
              
              
              %gain
                % rand seq
%               mParams = fitMultiVarParams(xBin([1,2],:)',yBin([1,2],:)',0,starterParams);
%               outNLGain = multiGainNL(mParams,xBin([1,2],:)');
%               
              
                % rand static
              
              mParams = fitMultiVarParams(xBin([1,3],:)',yBin([1,3],:)',0,starterParams);
              outNLGain = multiGainNL(mParams,xBin([1,3],:)');
              
              staticErrorGain(d) = immse(yBin(3,:)',outNLGain(:,2));
              MSErandGain2(d) = immse(yBin(1,:)',outNLGain(:,1));
              MSErandGain(d) = MSErandGain2(d)/staticErrorGain(d);
              
                % seq static
              mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',0,starterParams);
              outNLGain = multiGainNL(mParams,xBin([2,3],:)');
              

              staticErrorGain(d) = immse(yBin(3,:)',outNLGain(:,2));
              MSEseqGain2(d) = immse(yBin(2,:)',outNLGain(:,1));
              MSEseqGain(d) = MSEseqGain2(d)/staticErrorGain(d);
              
%               MSEdecreaseGain(d) = (MSErandGain-MSErandHZ);
%               MSEdecrease(d) = (MSEseqGain-MSEseqHZ);
%               
              
              
              %%%%%BOTH%%%%%
%                 % rand seq
%               mParams = fitMultiVarParams(xBin([1,2],:),yBin([1,2],:),2,starterParams);
%               outNLHZGain = multiVarNL(mParams,xBin([1,2],:));

                % rand static
              mParams = fitMultiVarParams(xBin([1,3],:)',yBin([1,3],:)',2,starterParams);
              outNLBoth = multiVarNL(mParams,xBin([1,3],:)');
              

              
%               RandHZ1(d) = mParams(5)-mParams(3);
%               RandGain1(d) = mParams(4)-mParams(2);
              tester(d,:)=mParams;
              
              RANDpChangeHZ(d)  = (mParams(5)-mParams(3))/mParams(5);
              RANDpChangeGain(d) = (mParams(4)-mParams(2))/(-1*mParams(4));
              %Negative denominator should put gain and hz shift into
              %register beause the horizontal shift numbers are always
              %negative for some reason (but leftward shift is
              %*positive*)
              
                % seq static
              mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',2,starterParams);
              outNLBoth = multiVarNL(mParams,xBin([2,3],:)');
              
              

              
              tester2(d,:)=mParams;
            
              SEQpChangeHZ(d)  = (mParams(5)-mParams(3))/mParams(5);
              SEQpChangeGain(d) = (mParams(4)-mParams(2))/(-1*mParams(4));
              %low gain = low parameter, but not dividing by a negative num
              
%               SeqHZ1(d) = mParams(5)-mParams(3);
%               SeqGain1(d) = mParams(4)-mParams(2);
              %percent change final-initial/initial
              

              
        end
   

        seqChange=[];
randChange=[];
seqChange = 100*(multiCellSeq-multiCellStatic)./multiCellStatic;
randChange = 100*(multiCellRand-multiCellStatic)./multiCellStatic;

seqX = ones(1,length(multiCellSeq));
figure
% seqX=seqX+1;
changeY = [seqChange' randChange']
changeX = [seqX' seqX'+.2]
plot(changeX, changeY,'.')
axis([.5 1.7 -20 80])
line([seqX(1)' seqX(1)'+.2],[seqChange(1:end)' randChange(1:end)'])
line([seqX(1)' seqX(1)'+.2],[mean(seqChange) mean(randChange)],'LineWidth',3)
% makeAxisStruct(gca,strtrim(['changeDataMikeLINES' 'offSmooth']))
        
%         multiCellSeq(isnan(multiCellSeq))=[];
%         multiCellRand(isnan(multiCellRand))=[];
%         multiCellStatic(isnan(multiCellStatic))=[];
       
        cellType = fileIndex(:,1:10,1);
        
%         RandHZ = RandHZ1(1:d)';
%         RandGain = RandGain1(1:d)';
%         SeqHZ = SeqHZ1(1:d)';
%         SeqGain = SeqGain1(1:d)';
        
        RandHZOut=RANDpChangeHZ(1:d)'*100;
        RandGainOut = RANDpChangeGain(1:d)'*100;
        SeqHZOut = SEQpChangeHZ(1:d)'*100;
        SeqGainOut = SEQpChangeGain(1:d)'*100;
%         meanDecrease = mean(MSEdecrease);
%         errorDecrease = sem(MSEdecrease);
        
        
%         save(strcat('x',cellType,'_MN'),'RandHZOut','RandGainOut','SeqHZOut','SeqGainOut','meanDecrease','errorDecrease')
%         save(strcat(cellType,'_MN'),'RandHZ','SeqHZ','RandGain','SeqGain','meanDecrease','errorDecrease')
% FOR error graphs
%     figure
%     XforGraph = ones(1,length(RandHZOut));
%     changeY = [MSEseqHZ' MSEseqGain'];
%     changeX = [XforGraph' XforGraph'+.2];
%     plot(changeX, changeY,'.')
%     axis([.5 1.7 -50 50])
%     line([XforGraph(1)' XforGraph(1)'+.2],[MSEseqHZ(1:end)' MSEseqGain(1:end)'])
%     line([XforGraph(1)' XforGraph(1)'+.2],[mean(MSEseqHZ) mean(MSEseqGain)],'LineWidth',3)
%     title('HZ vs Gain')
% 
%     figure
%     XforGraph = ones(1,length(RandHZOut));
%     changeY = [MSEseqHZ' MSErandHZ'];
%     changeX = [XforGraph' XforGraph'+.2];
%     plot(changeX, changeY,'.')
%     axis([.5 1.7 -50 50])
%     line([XforGraph(1)' XforGraph(1)'+.2],[MSEseqHZ(1:end)' MSErandHZ(1:end)'])
%     line([XforGraph(1)' XforGraph(1)'+.2],[mean(MSEseqHZ) mean(MSErandHZ)],'LineWidth',3)
%     title('HZ')
%     % 
%     %         
%     figure
%     XforGraph = ones(1,length(RandHZOut));
%     changeY = [MSEseqGain' MSErandGain'];
%     changeX = [XforGraph' XforGraph'+.2];
%     plot(changeX, changeY,'.')
%     axis([.5 1.7 -30 30])
%     line([XforGraph(1)' XforGraph(1)'+.2],[MSEseqGain(1:end)' MSErandGain(1:end)'])
%     line([XforGraph(1)' XforGraph(1)'+.2],[mean(MSEseqGain) mean(MSErandGain)],'LineWidth',3)
%     title('Gain')

%%%%%%%%%%%%%%%%%%% some bullshit unity line

meanErrorHZ = mean(MSEseqHZ);
semErrorHZ=sem(MSEseqHZ);

meanErrorGain = mean(MSEseqGain);
semErrorGain = sem(MSEseqGain);


figure
changeY = [1-(MSEseqHZ/max([MSEseqHZ MSEseqGain]))'];
changeX = [1-(MSEseqGain/max([MSEseqHZ MSEseqGain]))'];
plot(changeX, changeY,'.')
line([0 1],[0 1],'Color','k')
title(cellType)
hold on
changeY = [1-MSErandHZ/max([MSErandHZ MSErandGain])' ];
changeX = [1-MSErandGain/max([MSErandHZ MSErandGain])'];
plot(changeX, changeY,'.')

figure
changeY = [MSEseqHZ'];
changeX = [MSEseqGain'];
plot(changeX, changeY,'.')
hold on
changeY = [MSErandHZ'];
changeX = [MSErandGain'];
plot(changeX, changeY,'.')
line([0 max([MSErandGain MSErandHZ MSEseqHZ MSEseqGain])],[0 max([MSErandGain MSErandHZ MSEseqHZ MSEseqGain])],'Color','k')
title(strtrim([cellType ' Straight MSE']))
plot(meanErrorGain,meanErrorHZ)
line([(meanErrorGain-semErrorGain) (meanErrorGain+semErrorGain)], [meanErrorHZ meanErrorHZ])
line([(meanErrorGain) (meanErrorGain)], [meanErrorHZ-semErrorHZ meanErrorHZ+semErrorHZ])
% makeAxisStruct(gca,strtrim(['MikeUnity_' 'onSmoothNorm']))

meanMSESeqHZ(typeInd)=mean(MSEseqHZ);
semMSESeqHZ(typeInd)=sem(MSEseqHZ);
meanMSESeqGain(typeInd)=mean(MSEseqGain);
semMSESeqGain(typeInd)=sem(MSEseqGain);
meanMSErandHZ(typeInd)=mean(MSErandHZ);
meanMSErandGain(typeInd)=mean(MSErandGain);


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
