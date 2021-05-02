%% 
exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
protocolToAnalyze = 'Motion And Noise';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:3)=[]; %removes macOS hidden stuff
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
            
yBin(1,xBin(1,:)<0)=0;
yBin(2,xBin(2,:)<0)=0;
yBin(3,xBin(3,:)<0)=0;


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
              
           
            MSErandHZ = immse(yBin(1,:)',outNLHZ(:,1));

              
                % seq static
              mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',1,starterParams);
              outNLHZ = multiHZNL(mParams,xBin([2,3],:)');
              

              
            MSEseqHZ = immse(yBin(2,:)',outNLHZ(:,1));
              
              
              %gain
                % rand seq
%               mParams = fitMultiVarParams(xBin([1,2],:)',yBin([1,2],:)',0,starterParams);
%               outNLGain = multiGainNL(mParams,xBin([1,2],:)');
%               
              
                % rand static
              mParams = fitMultiVarParams(xBin([1,3],:)',yBin([1,3],:)',0,starterParams);
              outNLGain = multiGainNL(mParams,xBin([1,3],:)');
              
              
              MSErandGain = immse(yBin(1,:)',outNLGain(:,1));
              
                % seq static
              mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',0,starterParams);
              outNLGain = multiGainNL(mParams,xBin([2,3],:)');
              

              
              MSEseqGain = immse(yBin(2,:)',outNLGain(:,1));
              
              MSEdecreaseGain(d) = (MSErandGain-MSErandHZ);
              MSEdecrease(d) = (MSEseqGain-MSEseqHZ);
              
              
              
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
        cellType = fileIndex(:,1:10,1);
        
%         RandHZ = RandHZ1(1:d)';
%         RandGain = RandGain1(1:d)';
%         SeqHZ = SeqHZ1(1:d)';
%         SeqGain = SeqGain1(1:d)';
        
        RandHZOut=RANDpChangeHZ(1:d)'*100;
        RandGainOut = RANDpChangeGain(1:d)'*100;
        SeqHZOut = SEQpChangeHZ(1:d)'*100;
        SeqGainOut = SEQpChangeGain(1:d)'*100;
        meanDecrease = mean(MSEdecrease);
        errorDecrease = sem(MSEdecrease);
        
        
        save(strcat('x',cellType,'_MN'),'RandHZOut','RandGainOut','SeqHZOut','SeqGainOut','meanDecrease','errorDecrease')
%         save(strcat(cellType,'_MN'),'RandHZ','SeqHZ','RandGain','SeqGain','meanDecrease','errorDecrease')
   
% X=[];
% Y=[];
% for g= 1:3
% nlParams(g,:) = models.ln.fitNonlinearityParams(X(g,:), Y(g,:));
% plot(X(g,:),outputNonlinearity(nlParams(g,:),X(g,:)))
% hold on
% end
% X=X'
% Y=Y'

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
