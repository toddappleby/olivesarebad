 %% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = getDirectory()
% exportFolder = 'C:\Users\todda\Documents\Primate Data\ClarinetExports\';
protocolToAnalyze = 'Motion And Noise';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[]; %removes macOS hidden stuff .. sometimes 2 sometimes 3???
% folderSet(:,:,1:3)=[]; 
%%
 timings = [250,10000,250]; %same for these
 nonlinearityBins = 100;
kmeansGroups=[]; 
cellnameList = [];
errorGroups = [];


recordingType = 'extracellular';
binRate = 1e3;
spikeMeans = cell(13,size(folderSet,3));
% 
    allCellNames = [];   
    allCellTypes = [];
    allSeqMeans = [];
    allRandMeans = [];
    allIndex = [];
%   1:size(folderSet,3)
count=0;
createSharableFormat = 1;
% includeTypes =[13 12 8 7];
includeTypes = 13;
%
for c = includeTypes
 
    cd(strcat(typesFolder,folderSet(:,:,c))) 
    
    cellType = strtrim(folderSet(:,:,c)); 

    
    typeDir = dir; 
    fileSet = char(string({typeDir.name})); 
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    cellTypeData = [];
    seqSeedMean = [];
    randSeedMean = [];
    staticSeedMean = [];
    seqMean = [];
    randMean = [];
    staticMean = [];
    seqMeanFirst5 = [];
    randMeanFirst5 = [];
    staticMeanFirst5 = [];
    seqMeanLast5 = [];
    randMeanLast5 = [];
    staticMeanLast5 = [];
    MSErandGain=[];
    MSEseqGain=[];
%     RANDpChangeHZ=[];
%     RANDpChangeGain=[];
    MSEseqHZ =[];
    MSErandHZ=[];
    cellList = char;
    cellTypeList=string;
pcChangeHzSequential = [];
pcChangeHzRandom = [];
pcChangeGainSequential = [];
pcChangeGainRandom = [];
hzRSE = [];
gainRSE = [];

        for d = 1:size(fileIndex,3)
%           for d = 22:22
           outI = []; 
           
         
           
           
         yearRecorded = fileIndex(:,1:4,d);
         cellName = strtrim(fileIndex(:,:,d));
         cellList(d,:) = cellName;
         cellTypeList(d,:)=string(cellType);
%          if strcmp(yearRecorded,'2019')
%              continue
%          end
           
         load(strtrim(fileIndex(:,:,d)))
         %Where params are most appropriate?  prob the ones I ran the most:
         clear indLength
         
         
         
            for e = 1:size(indexHolder,2)
                indLength(e) = length(indexHolder{2,e});
            end 
            
            maxI = maxk(indLength,3); %top 3 conditions -- must be best options for seq rand & static
            logI = ismember(string(indLength),string(maxI));
            outI = find(logI==1); %index holder indices of top options
            matrixSplit=[];
            logicalMatrix=[];
%             matrixSplit = reshape([indexHolder{1,:}],[6,length(indexHolder)])';
            matrixSplit = vertcat(indexHolder{1,:});
            
           if strcmp(cellType,'ParasolOFF') 
             switch cellName
                 case '20200108Bc1_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "250", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                 case '20200504Bc1_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "300", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                 case '20210128Bc2_MN.mat'
                    logicalMatrix = contains(matrixSplit,["binary", "Control",  "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==4)];
                 case '20210507Bc4_MN.mat'
                    logicalMatrix = contains(matrixSplit,["binary", "Control",  "90", "stationary","sequential","random"])
                    outI = [find(sum(logicalMatrix,2)==4)];
                 case '20190923Bc2_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "300", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                 case '20191003Bc7_MN.mat'
                     logicalMatrix = contains(matrixSplit,["gaussian", "Control", "325", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                 case '20191125Bc4_MN.mat'
                      logicalMatrix = contains(matrixSplit,["gaussian", "Control","300", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                 case '20191126Bc3_MN.mat'
                      logicalMatrix = contains(matrixSplit,["gaussian", "Control", "350", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                 case '20220406Ac1_MN.mat'
                      logicalMatrix = contains(matrixSplit,["gaussian", "Control", "350", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                    
%                  case '20191105Bc3_MN.mat'
%                      logicalMatrix = contains(matrixSplit,["gaussian", "Control", "350", "90", "stationary","sequential","random"]);
%                     outI = [find(sum(logicalMatrix,2)==5)]
                    
             end
           end
           
           if strcmp(cellType,'SmoothON') 
             switch cellName
                 case '20200615Bc1_MN.mat' 
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "375", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                 case '20200630Ac1_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "motion and noise", "180", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==4)];
                 case '20210122Ac1_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "450", "180", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];      
                 case '20230320Ac5_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "375", "180", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];    
                 case '20230725Ac1_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "275", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];    
             end
           end
           
           if strcmp(cellType,'Lbs') 
             switch cellName
                 case '20220406Ac1_MN.mat'
                    logicalMatrix = contains(matrixSplit,["binary", "AltCenter", "450", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];    
             end
           end
           
           if strcmp(cellType,'BroadThorny') 
             switch cellName
                 case '20211102Gc8_MN.mat'
                    logicalMatrix = contains(matrixSplit,["gaussian", "Control", "525", "90", "stationary","sequential","random"]);
                    outI = [find(sum(logicalMatrix,2)==5)];
                    
             end
           end
            
            if strcmp(strtrim(fileIndex(:,:,d)),'20210507Bc4_MN.mat')
                outI = [1 2 3];
            end
            
            
            %which indices correspond to which condition?
            clear holdLabel
            for ee = 1:length(outI)
                holdLabel(ee) = string(indexHolder{1,outI(ee)}(size(matrixSplit,2))); %
            end
            
            
            seqIndex = indexHolder{2,outI(strcmp(holdLabel,"sequential"))};
            randIndex = indexHolder{2,outI(strcmp(holdLabel,"random"))};
            staticIndex = indexHolder{2,outI(strcmp(holdLabel,"stationary"))};
            
            
            
            prePts = 250 * 1e-3 * metaData.sampleRate;
            stimPts = 10000 * 1e-3 * metaData.sampleRate;
            tailPts = 250 * 1e-3 * metaData.sampleRate;
           lfilter=[];
           % going forward: grab all frame sequences at once.
           % use logical indexing to separate conditions out for linear
           % filter and NLs
            for f = 1:3
                indices = indexHolder{2,outI(f)};
                frameDwell = ones(length(indices),1)*double(indexHolder{1,outI(f)}(3));
                
                noiseVars = struct();
                noiseVars.contrast = 0.3333;
                noiseVars.type = indexHolder{1,outI(f)}(1);

                frames = manookinlab.ovation.getFrameTimesFromMonitor(frameTimings(indices+size(spikingData,1),:), metaData.sampleRate, binRate);
          if strcmp(cellType,'BroadThorny') && strcmp(cellName,'20230613Bc5_MN.mat') || strcmp(cellName,'20200713Bc3_MN.mat')
              cellName
                frameSeq = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),binRate,frames,4,seed(indices),frameDwell);
          else
                frameSeq = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),binRate,frames,2,seed(indices),frameDwell);
          end
                response = zeros(size(indices,1),metaData.stimTime+metaData.tailTime);
                for ff = 1:length(indices)
                    response(ff,:) = binSpikeCount(spikingData(indices(ff),prePts+1:end), binRate, metaData.sampleRate);
                    response(ff,:) = psth(response(ff,:),6+2/3,binRate,1);
                end
                
                
             stimulus = frameSeq(seed(indices) ~= 1, metaData.preTime+1:end);
             responses = response(seed(indices) ~= 1, :);
                
            lfilter(f,:) = getLinearFilter(stimulus, responses, ...
                'analysisType', 'revcorr', ...
                'fourierCorrection',false, ...
                'binRate', binRate, ...
                'filterTime', 0.5, ...
                'frameRate', metaData.frameRate);
            lfilter(f,:) = lfilter(f,:) / norm(lfilter(f,:));
            end
           allFilters = lfilter;
lfilter = mean(lfilter,1); %without mean filter, NLs sometimes look different
% figure(5)
% plot(lfilter(1:500),'m','LineWidth',2)
% hold on
% plot(allFilters(1,1:500),'--b')
% plot(allFilters(2,1:500),'--r')
% plot(allFilters(3,1:500),'--k')
% legend('AverageFilter','RandomFilter', 'SequentialFilter','StationaryFilter')
% 
% title('sta')
% xlabel('ms')
% ylabel('contrast weight')

            stimIndex = metaData.preTime + (1 : metaData.stimTime);
%             stimIndex = preTime + (1 : stimTime);
%            doing it twice because no time to rewrite simple thing?
            for g = 1:3
                
                indices = indexHolder{2,outI(g)};
                frameDwell = ones(length(indices),1)*double(indexHolder{1,outI(g)}(3));
                
                noiseVars = struct();
                noiseVars.contrast = 0.3333;
                noiseVars.type = indexHolder{1,outI(g)}(1);
                
                
                frames = manookinlab.ovation.getFrameTimesFromMonitor(frameTimings(indices+size(spikingData,1),:), metaData.sampleRate, binRate);
          if strcmp(cellType,'BroadThorny') && strcmp(cellName,'20230613Bc5_MN.mat') || strcmp(cellName,'20200713Bc3_MN.mat')
                frameSeq = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),binRate,frames,4,seed(indices),frameDwell);
          else
                frameSeq = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),binRate,frames,2,seed(indices),frameDwell);
          end                
                response = zeros(size(indices,1),metaData.stimTime+metaData.tailTime);
                for gg = 1:length(indices)
                    response(gg,:) = binSpikeCount(spikingData(indices(gg),prePts+1:end), binRate, metaData.sampleRate);
                    response(gg,:) = psth(response(gg,:),6+2/3,binRate,1);
                end
                 
%                 S = frameSeq;
%                 R = response;
                S = frameSeq(seed(indices) ~= 1, metaData.preTime+1:end);
                R = response(seed(indices) ~= 1, :);
                lf = lfilter(1 : size(S,2));
                P = zeros(size(R));
                for p = 1 : size(S,1)
                    tmp = ifft( fft(S(p,:)) .* fft(lf) );
                    P(p,:) = real(tmp);
                    P(p,:) = P(p,:)./std(P(p,:));
                end
                [xBin(g,:),yBin(g,:)] = binNonlinearity(P(:,stimIndex(:)),R(:,stimIndex(:)),nonlinearityBins);
                
classLabels = {'random', 'motion', 'stationary'};
responseSet{1,g} = classLabels{g};                 
responseSet{2,g} = response;      
%                 nlParams(g,:) = models.ln.fitNonlinearityParams(xBin(g,:), yBin(g,:));
%                 prediction(indices,:) = P;
            end
         
   %To create common save w/ Mike's format
   if createSharableFormat == true

    analysisParams = struct;
    stimulusParams = struct;

    analysisParams.nonlinearityBins = nonlinearityBins;
    analysisParams.xBin = xBin;
    analysisParams.yBin = yBin;
    analysisParams.lfilter = lf;
    analysisParams.lfilters = allFilters;
    analysisParams.classLabels = {'random', 'motion', 'stationary'};

    stimulusParams.seeds = seed;
    stimulusParams.backgroundClasses = splitCell{2,size(splitCell,2)};
    stimulusParams.response = responseSet;
    stimulusParams.stimulus = frameSeq;
    stimulusParams.sampleRate = metaData.sampleRate;
    stimulusParams.stimTime = metaData.stimTime;
    stimulusParams.preTime = metaData.preTime;
    stimulusParams.frameRate = metaData.frameRate;
    stimulusParams.cellType = cellType;
    stimulusParams.cellName = cellName(1:11);
    stimulusParams.frameTimings = frames;
    stimulusParams.tailTime = metaData.tailTime;
    stimulusParams.prediction = P;
    

    save(strcat(cellType,'_',cellName(1:11),'_MN_Sharable.mat'),'analysisParams','stimulusParams')
   end

%             mean of stimulus from each cycle
%                 figure(45)    
%       plot(lfilter(1:1000))      
%       figure(46)
%       plot(xBin',yBin')
      
% LF and NL subplot for each cell of each type
figure(100)
subplot(ceil(size(fileIndex,3)/2),2,d)

plot(xBin(1,:),yBin(1,:),'b');hold on;plot(xBin(2,:),yBin(2,:),'r');plot(xBin(3,:),yBin(3,:),'k')
title(cellName(1:11))

figure(99)
subplot(ceil(size(fileIndex,3)/2),2,d)

plot(lf(1:650))
title(cellName(1:11))
            
             
            
% spike count starts here  %%%%%%%%%%%%%%%%%%%%%%


%repeated seeds

seedInd = find(seed==1);   
I=[];
II=[];
III=[];
[I,J]=find(seqIndex==seedInd); %i dunno how to do logical indexing by row/column rather than the linearly (which creates out of array bounds error)
seqSeedMean(d) = mean(mean(sum(spikingData(I,:),2),2));

[II,J]=find(randIndex==seedInd);
randSeedMean(d) = mean(mean(sum(spikingData(II,:),2),2));

[III,J]=find(staticIndex==seedInd);
staticSeedMean(d) = mean(mean(sum(spikingData(III,:),2),2));



%all seeds
seqMean(d) = mean(mean(sum(spikingData(seqIndex(seed(seqIndex)~=1),5000:end),2),2));
randMean(d) = mean(mean(sum(spikingData(randIndex(seed(randIndex)~=1),5000:end),2),2));
staticMean(d) = mean(mean(sum(spikingData(staticIndex(seed(staticIndex)~=1),5000:end),2),2));

seqMeanInitial(d) = mean(mean(sum(spikingData(seqIndex,5000:10000),2),2));
randMeanInitial(d) = mean(mean(sum(spikingData(randIndex,5000:10000),2),2));
staticMeanInitial(d) = mean(mean(sum(spikingData(staticIndex,5000:10000),2),2));


seqMeanFirst5(d) = mean(mean(sum(spikingData(seqIndex,1000:50000),2),2));
randMeanFirst5(d) = mean(mean(sum(spikingData(randIndex,1000:50000),2),2));
staticMeanFirst5(d) = mean(mean(sum(spikingData(staticIndex,1000:50000),2),2));

seqMeanLast5(d) = mean(mean(sum(spikingData(seqIndex,50000:100000),2),2));
randMeanLast5(d) = mean(mean(sum(spikingData(randIndex,50000:100000),2),2));
staticMeanLast5(d) = mean(mean(sum(spikingData(staticIndex,50000:100000),2),2));
            
            
            
% fits here %%%%%%%%%%%%%%%%%%%%            
                %go with 1 being random, 2 sequential, 3 static
              starterParams = [.2, 0];
              %horizontal  
                % All
              %each half modeled
              if strcmp(cellType,'BroadThorny')
                 mParams = fitMultiVarParams(xBin([1,2,3],51:100),yBin([1,2,3],51:100),1,starterParams);
                 outNLHZOn = multiHZNL(mParams,xBin([1,2,3],51:100));

                 figure;colororder(['b';'r';'k']);plot(xBin(:,51:100)',outNLHZOn','--');hold on;plot(xBin(:,51:100)',yBin(:,51:100)')
                 

                        mParams2 = fitMultiVarParams(xBin([1,2,3],1:50),yBin([1,2,3],1:50),1,starterParams);
                 outNLHZOff = multiHZNL(mParams2,xBin([1,2,3],1:50));
                 
                 

                 colororder(['b';'r';'k']);plot(xBin(:,1:50)',outNLHZOff','--');hold on;plot(xBin(:,1:50)',yBin(:,1:50)')
                 title('Broad Thorny NLs w/ Separate On/Off Model (XshiftOnly)')
                 xlabel('Spike Rate (Hz)')
                 ylabel('Input')
                
                 mParams = mParams2;
                 
              else

              mParams = fitMultiVarParams(xBin([1,2,3],:),yBin([1,2,3],:),1,starterParams);
              outNLHZ = multiHZNL(mParams,xBin([1,2,3],:));
              
              end
              hzRSE(d,:) = [mParams(3),mParams(4),mParams(5)];

%                  mParams = fitMultiVarParams(xBin([1,2,3],51:100),yBin([1,2,3],51:100),1,starterParams);
%                  outNLGainOn = multiGainNL(mParams,xBin([1,2,3],51:100)');
% 
%                  figure;colororder(['b';'r';'k']);plot(xBin(:,51:100)',outNLGainOn,'--');hold on;plot(xBin(:,51:100)',yBin(:,51:100)')
%                  
% 
%                         mParams2 = fitMultiVarParams(xBin([1,2,3],1:50),yBin([1,2,3],1:50),1,starterParams);
%                  outNLGainOff = multiGainNL(mParams2,xBin([1,2,3],1:50)');
% 
%                  colororder(['b';'r';'k']);plot(xBin(:,1:50)',outNLGainOff,'--');hold on;plot(xBin(:,1:50)',yBin(:,1:50)')
%                  title('Broad Thorny NLs w/ Separate On/Off Model (Gain shift Only)')
%                  xlabel('Spike Rate (Hz)')
%                  ylabel('Input')
%                 pause



%try inverting the x shift params only (model both sides..)
%               if strcmp(cellType,'BroadThorny')
%                  mParams2 = fitMultiVarParams(xBin([1,2,3],1:50),yBin([1,2,3],1:50),1,starterParams);
%                  outNLHZOff = multiHZNL(mParams2,xBin([1,2,3],1:50));
% 
%                  colororder(['b';'r';'k']);plot(xBin(:,1:50)',outNLHZOff','--');hold on;plot(xBin(:,1:50)',yBin(:,1:50)')
% 
% 
%                  mParams = fitMultiVarParams(xBin([1,2,3],51:100),yBin([1,2,3],51:100),1,starterParams);
%                  outNLHZOn = multiHZNL([mParams(1) mParams(2) mParams(3:5)],xBin([1,2,3],51:100));
% 
%                  figure;colororder(['b';'r';'k']);plot(xBin(:,51:100)',-outNLHZOff','--');hold on;plot(xBin(:,51:100)',yBin(:,51:100)')
%                  
% 
%                     
% 
%                  colororder(['b';'r';'k']);plot(xBin(:,1:50)',outNLHZOff','--');hold on;plot(xBin(:,1:50)',yBin(:,1:50)')
%                  title('Broad Thorny NLs w/ Separate On/Off Model (XshiftOnly)')
%                  xlabel('Spike Rate (Hz)')
%                  ylabel('Input')
%                 pause
              %whole model then On half flipped params
%               if strcmp(cellType,'BroadThorny')
%                  mParams = fitMultiVarParams(xBin([1,2,3],1:100),yBin([1,2,3],1:100),1,starterParams);
%                  outNLHZOn = multiHZNL([-mParams(1:2),-mParams(3:5)],xBin([1,2,3],51:100));
%                  outNLHZOff = multiHZNL(mParams,xBin([1,2,3],1:50));
% 
%                  outNLHZ = [outNLHZOff,outNLHZOn];
% 
%                  figure;colororder(['b';'r';'k']);plot(xBin',outNLHZ','--');hold on;plot(xBin',yBin')
%                  title('Broad Thorny NLs w/ Flipped On/Off Model (XshiftOnly)')
%                  xlabel('Spike Rate (Hz)')
%                  ylabel('Input')                 
% pause


     

%               figure(20)
%               colororder(['b';'r';'k']) %matches above order
%               plot(xBin',outNLHZ','--')
%               hold on
%               plot(xBin',yBin')
%               

              
              
%               figure(19)
%               plot(xBin(1,:),outNLHZ(1,:),'b')
%               hold on
%               plot(xBin(2,:),outNLHZ(2,:),'r')
%               plot(xBin(3,:),outNLHZ(3,:),'k')
%          
%               plot(xBin(1,:),yBin(1,:),'b--')
%               plot(xBin(2,:),yBin(2,:),'r--')
%               plot(xBin(3,:),yBin(3,:),'k--')
%               legend('Random-Model','Sequential-Model','Random-Data','Sequential-Data')
%               title('Raw and Modeled NLs - XShift')
%               xlabel('input')
%               ylabel('spike rate (Hz)')
%               
              

              
% %                rand static
%               mParams = fitMultiVarParams(xBin([1,3],:),yBin([1,3],:),1,starterParams);
%               outNLHZ = multiHZNL(mParams,xBin([1,3],:));
%               
%               
%             staticErrorHZ1(d) = immse(yBin(3,:),outNLHZ(2,:));  
%             MSErandHZ2(d) = immse(yBin(1,:),outNLHZ(1,:));
%             MSErandHZ(d) = MSErandHZ2(d)/staticErrorHZ1(d);
%             RSQrandHZ(d) = R2(yBin(1,:)',outNLHZ(1,:));
%             
%               
%                 % seq static
%               mParams = fitMultiVarParams(xBin([2,3],:),yBin([2,3],:),1,starterParams);
%               outNLHZ = multiHZNL(mParams,xBin([2,3],:));
%             
%             staticErrorHZ2(d) = immse(yBin(3,:),outNLHZ(2,:));  
%             MSEseqHZ2(d) = immse(yBin(2,:),outNLHZ(1,:));
%             MSEseqHZ(d) = MSEseqHZ2(d)/staticErrorHZ2(d);
%             RSQseqHZ(d) = R2(yBin(1,:)',outNLHZ(1,:));
%               
%             
% [rsquaredOut,rootMSE]=rsquare(yBin(2,:)',outNLHZ(:,1))
% 
% [rsquaredOut,rootMSE]=rsquare(binnedY(:,2),outNL(:,2))
% 
% MSE1 = immse(binnedY(:,1),outNL(:,1))
% MSE2 = immse(binnedY(:,2),outNL(:,2))
            
              
              %gain
                % gain all
              mParams = fitMultiVarParams(xBin,yBin,0,starterParams);
              outNLGain = multiGainNL(mParams,xBin);
              gainRSE(d,:) = [mParams(3), mParams(4), mParams(5)];
              
%               figure(21)
%               colororder(['b';'r';'k']) %matches above order
%               plot(xBin',outNLGain,'--')
%               hold on
%               plot(xBin',yBin')
%               legend('Random','Seq','Static')
              
%               figure(18)
%            
%               plot(xBin(1,:),outNLGain(:,1),'b')
%               hold on
%               plot(xBin(2,:),outNLGain(:,2),'r')
%               plot(xBin(3,:),outNLGain(:,2),'k')
%               plot(xBin(1,:),yBin(1,:),'b--')
%               plot(xBin(2,:),yBin(2,:),'r--')
%               plot(xBin(3,:),yBin(3,:),'k--')
%               legend('Random-Model','Sequential-Model','Random-Data','Sequential-Data')
%               title('Raw and Modeled NLs - Gain Shift')
%               xlabel('input')
%               ylabel('spike rate (Hz)')
%               
%              
%               
%               
%                 %rand static
%               mParams = fitMultiVarParams(xBin([1,3],:)',yBin([1,3],:)',0,starterParams);
%               outNLGain = multiGainNL(mParams,xBin([1,3],:)');
%               
%               staticErrorGain1(d) = immse(yBin(3,:)',outNLGain(:,2));
%               MSErandGain2(d) = immse(yBin(1,:)',outNLGain(:,1));
%               MSErandGain(d) = MSErandGain2(d)/staticErrorGain1(d);
%               RSQrandGain(d) = R2(yBin(1,:)',outNLGain(:,1));
%               
%              
%               
%                 % seq static
%               mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',0,starterParams);
%               outNLGain = multiGainNL(mParams,xBin([2,3],:)');
%               
%               staticErrorGain2(d) = immse(yBin(3,:)',outNLGain(:,2));
%               MSEseqGain2(d) = immse(yBin(2,:)',outNLGain(:,1));
%               MSEseqGain(d) = MSEseqGain2(d)/staticErrorGain2(d); 
%               RSQseqGain(d) = R2(yBin(2,:)',outNLGain(:,1));
%               
%               MSEdecrease(d) = (MSEseqGain-MSEseqHZ)/MSEseqGain;
              
              
              
              %both 
%                 % rand seq
%               mParams = fitMultiVarParams(xBin([1,2],:),yBin([1,2],:),2,starterParams);
%               outNLHZGain = multiVarNL(mParams,xBin([1,2],:));

                % rand static
%               mParams =
%               fitMultiVarParams(xBin([1,3],:)',yBin([1,3],:)',2,starterParams); 
%               outNLBoth = multiVarNL(mParams,xBin([1,3],:)');
%               
%             RANDpChangeHZ(d)  = (mParams(5)-mParams(3))/mParams(5);
%               RANDpChangeGain(d) = (mParams(4)-mParams(2))/(-1*mParams(4));
%               
%               
%                 % seq static
%               mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',2,starterParams);
%               outNLBoth = multiVarNL(mParams,xBin([2,3],:)');
%               
%             
%              SEQpChangeHZ(d)  = (mParams(5)-mParams(3))/mParams(5);
%               SEQpChangeGain(d) = (mParams(4)-mParams(2))/(-1*mParams(4));
              %percent change final-initial/initial
  end              
            allTimes = [];
            spikeMatrix=[];
            

            for g = 1:size(spikingData,1)
                 spikes = binSpikeCount(spikingData(g,:), binRate, 10000);
                 spikes = spikes';
                 spikeMatrix(g,:) = spikes;
            end
            
            for t = staticIndex' 
                spikeTimes = find(spikeMatrix(t,:)>0);
                isiTmp = diff(spikeTimes);
                allTimes = [allTimes; isiTmp(:)];
            end

            [y,x]= getSpikeAutocorrelation(allTimes); 
            rowNumber = ceil(size(fileIndex,3)/2);
            figure(50+c)
            hold on
%             subplot(2,rowNumber,d)
            sgtitle(cellType)
            plot(x,y)
            
            %          
           
        
            
      
        
        % just look at the curves 
        

        
        
seqMean(seqMean==0)=[];
randMean(randMean==0)=[];
staticMean(staticMean==0)=[];
        
%         
% colors = [1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;.5 .2 .3;.7 .7 .7;.1 .2 .3;.4 .1 .8;.2 .3 .6;0 .1 .9;.4 .9 .1;.1 .1 .1;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;.5 .2 .3;.7 .7 .7;.1 .2 .3;.4 .1 .8;.4 .1 .9];
seqX = ones(1,length(seqMean));
randX = 2*ones(1,length(randMean));
staticX = 3*ones(1,length(staticMean));
% 
% scatter([seqX randX staticX],[seqMean randMean staticMean])
% % 
% % figure
% % plot(seqMean)
% % hold on
% % plot(randMean)
% % plot(staticMean)
% 
%  
seqChange = 100*(seqMean-staticMean)./staticMean;
randChange = 100*(randMean-staticMean)./staticMean;
meanSeqChange(c)=mean(seqChange);
meanRandChange(c)=mean(randChange);

% motionAdaptIndex = (seqMean-staticMean)./(seqMean+staticMean);

% 
% allSeqMeans = [allSeqMeans; seqChange'];
% allRandMeans = [allRandMeans; randChange'];
% allCellNames = [allCellNames; cellList];
% allCellNames2 = allCellNames(:,1:11,:);
% allCellTypes = [allCellTypes; cellTypeList];
% allIndex = [allIndex; motionAdaptIndex'];

% 
% 
% figure
% plot(seqChange)
% hold on
% plot(randChange)
% 
% % makeAxisStruct(gca,strtrim(['changeDataMike' folderSet(:,:,c)]))
% 
% figure
% scatter([1 2 3],[mean(seqMean) mean(randMean) mean(staticMean)],'SizeData',20,'LineWidth',4)
% 
% figure
% for j = 1:length(seqMean)
%     scatter([1 2 3],[seqMean(j) randMean(j) staticMean(j)],'MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:))
%     hold on
% end
% scatter([1 2 3],[mean(seqMean) mean(randMean) mean(staticMean)],'d','SizeData',20,'LineWidth',4)
% title(cellType)
% 
% 
figure
% seqX=seqX+1;    
changeY = [seqChange' randChange']
changeX = [seqX' seqX'+.2]
plot(changeX, changeY,'.')
axis([.5 1.7 -50 100])
line([seqX(1)' seqX(1)'+.2],[seqChange(1:end)' randChange(1:end)'])
line([seqX(1)' seqX(1)'+.2],[mean(seqChange) mean(randChange)],'LineWidth',3)
makeAxisStruct(gca,strtrim(['changeDataLINES' folderSet(:,:,c)]))

% figure(89+c)
% 
% maxMean = max([mean(seqMeanFirst5) mean(randMeanFirst5) mean(staticMeanFirst5) mean(seqMeanLast5) mean(randMeanLast5) mean(staticMeanLast5)]);


% scatter([1 2 3],[mean(seqMeanFirst5)/maxMean mean(randMeanFirst5)/maxMean mean(staticMeanFirst5)/maxMean],'*','SizeData',20,'LineWidth',4)
% hold on
% scatter([4 5 6],[mean(seqMeanLast5)/maxMean mean(randMeanLast5)/maxMean mean(staticMeanLast5)/maxMean],'*','SizeData',20,'LineWidth',4)
% title(cellType)
% plot(sem(seqMeanFirst5))
% plot(sem(randMeanFirst5))
% plot(sem(staticMeanFirst5))
% plot(sem(seqMeanLast5))
% plot(sem(randMeanLast5))
% plot(sem(staticMeanLast5))
% makeAxisStruct(gca,strtrim(['earlyLate_' folderSet(:,:,c)]))
% 
% sem(seqMeanFirst5/max(seqMeanFirst5))
% sem(randMeanFirst5/max(randMeanFirst5))
% sem(staticMeanFirst5/max(staticMeanFirst5))
% 
% sem(seqMeanLast5/max(seqMeanLast5))
% sem(randMeanLast5/max(randMeanLast5))
% sem(staticMeanLast5/max(staticMeanLast5))

%  seqSeedMean = seqSeedMean(~isnan(seqSeedMean));
%  randSeedMean = randSeedMean(~isnan(randSeedMean));
%  staticSeedMean=staticSeedMean(~isnan(staticSeedMean));
%         
%  sSeqMeanCells = mean(seqSeedMean); 
%  sSeqSEMCells = sem(seqSeedMean);
%  SeqMeanCells = mean(seqMean);
%  SeqSEMCells = sem(seqMean);
%  
%  sRandMeanCells = mean(randSeedMean); 
%  sRandSEMCells= sem(randSeedMean);
%  RandMeanCells = mean(randMean);
%  RandSEMCells= sem(randMean);
%  
%  sStaticMeanCells= mean(staticSeedMean);
%  sStaticSEMCells = sem(staticSeedMean);
%  StaticMeanCells = mean(staticMean);
%  StaticSEMCells = sem(staticMean);
%  
% spikeMeans{1,c} = strtrim(folderSet(1,1:size(folderSet,2),c));
% spikeMeans{2,c} = sSeqMeanCells; 
% spikeMeans{3,c} = sSeqSEMCells;
% spikeMeans{8,c} = SeqMeanCells;
% spikeMeans{9,c} = SeqSEMCells;
% spikeMeans{4,c} = sRandMeanCells;
% spikeMeans{5,c} = sRandSEMCells;
% spikeMeans{10,c} = RandMeanCells;
% spikeMeans{11,c} = RandSEMCells;
% spikeMeans{6,c} = sStaticMeanCells;
% spikeMeans{7,c} = sStaticSEMCells;
% spikeMeans{12,c} =  StaticMeanCells;
% spikeMeans{13,c} = StaticSEMCells;



%         
%         
%         RandHZOut=RANDpChangeHZ(1:d)'*100;
%         RandGainOut = RANDpChangeGain(1:d)'*100;
%         SeqHZOut = SEQpChangeHZ(1:d)'*100;
%         SeqGainOut = SEQpChangeGain(1:d)'*100;
%         meanDecrease = mean(MSEdecrease)*100;
%         errorDecrease = sem(MSEdecrease)*100;
        
%         MSErandGain
%         MSEseqGain
%         
%  figure
% XforGraph = ones(1,length(RandHZOut));
% changeY = [SeqHZOut RandHZOut]
% changeX = [XforGraph' XforGraph'+.2]
% plot(changeX, changeY,'.')
% axis([.5 1.7 -50 50])
% line([XforGraph(1)' XforGraph(1)'+.2],[SeqHZOut(1:end) RandHZOut(1:end)])
% line([XforGraph(1)' XforGraph(1)'+.2],[mean(SeqHZOut) mean(RandHZOut)],'LineWidth',3)
% title('HZ')


%%%%%%%%%%%%%%%%%%% these are unity line plots 


% load('/Users/reals/Documents/PhD 2021/ClarinetExports/MeansandErrorsofError_Mike')

% meanErrorHZ = mean(MSEseqHZ);
% semErrorHZ=sem(MSEseqHZ);
% meanErrorGain = mean(MSEseqGain);
% semErrorGain = sem(MSEseqGain);
% 
% meanErrorHZ = mean([meanMSESeqHZ(c-4) meanErrorHZ]);
% semErrorHZ = sem([semMSESeqHZ(c-4) semErrorHZ]);
% meanErrorGain = mean([meanMSESeqGain(c-4) meanErrorGain]);
% semErrorGain = sem([semMSESeqGain(c-4) semErrorGain]);
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
% 
% 
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
% 
% if c == 4
% makeAxisStruct(gca,strtrim(['ToddUnity_' 'offParasolnorm']))
% elseif c == 5
% makeAxisStruct(gca,strtrim(['ToddUnity_' 'offSmoothnorm']))
% elseif c == 6   
% makeAxisStruct(gca,strtrim(['ToddUnity_' 'onParasolnorm']))
% elseif c == 7    
% makeAxisStruct(gca,strtrim(['ToddUnity_' 'onSmoothnorm'])) 
% end
%     
% 
% 
% figure
% changeY = [RSQseqHZ];
% changeX = [RSQseqGain];
% plot(changeX, changeY,'.')
% line([0 1],[0 1],'Color','k')
% title(strtrim([cellType ' Rsquare']))
count=count+1;
      figure(29) % seq static
       pcChangeHzSequential =  (diff(hzRSE(:,[2 3]),1,2)./sum(hzRSE(:,[2 3]),2))*100;
       plot(count,mean(pcChangeHzSequential),'r.','MarkerSize',30)
       hold on
       plot(count,pcChangeHzSequential,'ro', 'MarkerSize',5)


       
     seqHzError(count) = sem(pcChangeHzSequential);

       pcChangeHzRandom =  (diff(hzRSE(:,[1 3]),1,2)./sum(hzRSE(:,[1 3]),2))*100;
       plot(count+.2,mean(pcChangeHzRandom),'b.','MarkerSize',30)
       plot(count+.2,pcChangeHzRandom,'bo', 'MarkerSize',5)
       xlabel('cell type')
       ylabel('Percent Change Horizontal Shift Parameter')
       title('Change in X Shift Parameter in NL Models for Motion & Random Conditions')
       
    randHzError(count) = sem(pcChangeHzRandom);
    

    
       
      figure(31) %seq static
       pcChangeGainSequential =  (diff(gainRSE(:,[2 3]),1,2)./sum(gainRSE(:,[2 3]),2))*100;
       plot(count,mean(pcChangeGainSequential),'r.','MarkerSize',30)
       hold on
%        plot(count,pcChangeGainSequential,'ro','MarkerSize',5)
       
       
    seqGainError(count) = sem(pcChangeGainSequential);  
       
      %rand static
       pcChangeGainRandom =  (diff(gainRSE(:,[1 3]),1,2)./sum(gainRSE(:,[1 3]),2))*100;
       plot(count+.2,mean(pcChangeGainRandom),'b.','MarkerSize',30)
%        plot(count+.2,pcChangeHzRandom,'bo', 'MarkerSize',5)
       xlabel('cell type')
       ylabel('Percent Change Gain Shift Parameter from Static Condition')
       title('Change in Gain Parameter in NL Models for Motion & Random Conditions')

           
    randGainError(count) = sem(pcChangeGainRandom);
       
     
    



end
%ERROR Line graphs
%       figure(30)
%        plot(seqHzError)
%        hold on
%        plot(randHzError)
% 
%       figure(33)
%        plot(seqGainError)
%        hold on
%        plot(randGainError)

  
% 
%         
% figure
% XforGraph = ones(1,length(RandHZOut));
% changeY = [SeqGainOut RandGainOut]
% changeX = [XforGraph' XforGraph'+.2]
% plot(changeX, changeY,'.')
% axis([.5 1.7 -30 30])
% line([XforGraph(1)' XforGraph(1)'+.2],[SeqGainOut(1:end) RandGainOut(1:end)])
% line([XforGraph(1)' XforGraph(1)'+.2],[mean(SeqGainOut) mean(RandGainOut)],'LineWidth',3)
% title('Gain')
% makeAxisStruct(gca,strtrim(['changeDataLINES' folderSet(:,:,c)]))


%%%%%%%  simple error graphs
% figure
% XforGraph = ones(1,length(RandHZOut));
% changeY = [MSEseqHZ' MSEseqGain'];
% changeX = [XforGraph' XforGraph'+.2];
% plot(changeX, changeY,'.')
% axis([.5 1.7 -50 50])
% line([XforGraph(1)' XforGraph(1)'+.2],[MSEseqHZ(1:end)' MSEseqGain(1:end)'])
% line([XforGraph(1)' XforGraph(1)'+.2],[mean(MSEseqHZ) mean(MSEseqGain)],'LineWidth',3)
% title('HZ vs Gain')

% figure
% XforGraph = ones(1,length(RandHZOut));
% changeY = [MSEseqHZ' MSErandHZ'];
% changeX = [XforGraph' XforGraph'+.2];
% plot(changeX, changeY,'.')
% axis([.5 1.7 -50 50])
% line([XforGraph(1)' XforGraph(1)'+.2],[MSEseqHZ(1:end)' MSErandHZ(1:end)'])
% line([XforGraph(1)' XforGraph(1)'+.2],[mean(MSEseqHZ) mean(MSErandHZ)],'LineWidth',3)
% title('HZ')
% % 
% %         
% figure
% XforGraph = ones(1,length(RandHZOut));
% changeY = [MSEseqGain' MSErandGain'];
% changeX = [XforGraph' XforGraph'+.2];
% plot(changeX, changeY,'.')
% axis([.5 1.7 -30 30])
% line([XforGraph(1)' XforGraph(1)'+.2],[MSEseqGain(1:end)' MSErandGain(1:end)'])
% line([XforGraph(1)' XforGraph(1)'+.2],[mean(MSEseqGain) mean(MSErandGain)],'LineWidth',3)
% title('Gain')
        
%         save(strcat(cellType,'_MN'),'RandHZOut','RandGainOut','SeqHZOut','SeqGainOut','meanDecrease','errorDecrease')
        
        
%         spikeMeans

% meanMSESeqHZ(c-3)=mean(MSEseqHZ);
% semMSESeqHZ(c-3)=sem(MSEseqHZ);
% meanMSESeqGain(c-3)=mean(MSEseqGain);
% semMSESeqGain(c-3)=sem(MSEseqGain);
% meanMSErandHZ(c-3)=mean(MSErandHZ);
% meanMSErandGain(c-3)=mean(MSErandGain);



% end     
% seqChangeTodd = seqChange;
% randChangeTodd = randChange;
% xforplot = [1 2 3 4];
% % save('figure1Todd','seqChangeTodd','randChangeTodd')
% plot(xforplot,flip(meanMSESeqHZ))
% hold on
% errbar(xforplot,flip(meanMSESeqHZ),flip(semMSESeqHZ))
% 
% plot(xforplot,flip(-meanMSESeqGain))
% hold on
% errbar(xforplot,flip(-meanMSESeqGain),flip(semMSESeqGain))
% plot([1 2 3 4],(flip(meanMSESeqHZ-meanMSESeqGain)))

% % makeAxisStruct(gca,strtrim(['HZGAINcomparison']))

% save('MeansandErrorsofError_Todd','meanMSESeqHZ','semMSESeqHZ','meanMSESeqGain','semMSESeqGain')
% cd(exportFolder)
%  save('PercentChangeMotionandNoise','allIndex','allSeqMeans','allRandMeans','allCellNames2','allCellTypes')
%% test 
mean([meangainTodd;meanMSESeqGain])
mean([semGainTodd;semMSESeqGain])
mean([meanHZTodd;meanMSESeqHZ])
mean([semHZTodd;semMSESeqHZ])

mean([meanHZTodd;meanMSESeqHZ])-mean([meangainTodd;meanMSESeqGain])

% plot(seqX,seqMean,'.')

% plot(seqX+.2,randMean,'.')
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
%% graphing 2

%BT

figure
test = bar(categorical({'Repeated','Random'}),[spikeMeans{2,1} spikeMeans{4,1} spikeMeans{6,1}; spikeMeans{8,1},spikeMeans{10,1},spikeMeans{12,1}],'Facecolor','flat');

test(1).CData(1,:)=[1 0 0];
test(1).CData(2,:)=[1 0 0];

test(2).CData(1,:)=[0 0 1];
test(2).CData(2,:)=[0 0 1];

test(3).CData(1,:)=[0 0 0];
test(3).CData(2,:)=[0 0 0];
ylabel('Mean Spike Count')
set(gca,'box','off') 
set(gca,'FontSize',20)
legend('Motion','Random','Static','Location','NorthEast')
title('Broad Thorny')



% OFF P


figure
test = bar(categorical({'Repeated','Random'}),[spikeMeans{2,3} spikeMeans{4,3} spikeMeans{6,3}; spikeMeans{8,3},spikeMeans{10,3},spikeMeans{12,3}],'Facecolor','flat');

test(1).CData(1,:)=[1 0 0];
test(1).CData(2,:)=[1 0 0];

test(2).CData(1,:)=[0 0 1];
test(2).CData(2,:)=[0 0 1];

test(3).CData(1,:)=[0 0 0];
test(3).CData(2,:)=[0 0 0];
ylabel('Mean Spike Count')
set(gca,'box','off') 
set(gca,'FontSize',20)
legend('Motion','Random','Static','Location','NorthEast')
title('Off Parasol')

% OFF Smooth

figure
test = bar(categorical({'Repeated','Random'}),[spikeMeans{2,4},spikeMeans{4,4},spikeMeans{6,4}; spikeMeans{8,4},spikeMeans{10,4},spikeMeans{12,4}],'Facecolor','flat');

test(1).CData(1,:)=[1 0 0];
test(1).CData(2,:)=[1 0 0];

test(2).CData(1,:)=[0 0 1];
test(2).CData(2,:)=[0 0 1];

test(3).CData(1,:)=[0 0 0];
test(3).CData(2,:)=[0 0 0];
ylabel('Mean Spike Count')
set(gca,'box','off') 
set(gca,'FontSize',20)
legend('Motion','Random','Static','Location','NorthEast')
title('Off Smooth')

% ON Parasol

figure
test = bar(categorical({'Repeated','Random'}),[spikeMeans{2,5},spikeMeans{4,5},spikeMeans{6,5}; spikeMeans{8,5},spikeMeans{10,5},spikeMeans{12,5}],'Facecolor','flat');

test(1).CData(1,:)=[1 0 0];
test(1).CData(2,:)=[1 0 0];

test(2).CData(1,:)=[0 0 1];
test(2).CData(2,:)=[0 0 1];

test(3).CData(1,:)=[0 0 0];
test(3).CData(2,:)=[0 0 0];
ylabel('Mean Spike Count')
set(gca,'box','off') 
set(gca,'FontSize',20)
legend('Motion','Random','Static','Location','NorthEast')
title('ON Parasol')

%ON Smooth

figure
test = bar(categorical({'Repeated','Random'}),[spikeMeans{2,6},spikeMeans{4,6},spikeMeans{6,6}; spikeMeans{8,6},spikeMeans{10,6},spikeMeans{12,6}],'Facecolor','flat');

test(1).CData(1,:)=[1 0 0];
test(1).CData(2,:)=[1 0 0];

test(2).CData(1,:)=[0 0 1];
test(2).CData(2,:)=[0 0 1];

test(3).CData(1,:)=[0 0 0];
test(3).CData(2,:)=[0 0 0];
ylabel('Mean Spike Count')
set(gca,'box','off') 
set(gca,'FontSize',20)
legend('Motion','Random','Static','Location','NorthEast')
title('ON Smooth')



% 
% % figure
% % bar([spikeMeans{8,1},spikeMeans{10,1},spikeMeans{12,1}])
% 
% figure
% bar([spikeMeans{2,3},spikeMeans{4,3},spikeMeans{6,3}])
% figure
% bar([spikeMeans{8,3},spikeMeans{10,3},spikeMeans{12,3}])
% 
% figure
% bar([spikeMeans{2,4},spikeMeans{4,4},spikeMeans{6,4}])
% title('OFF Smooth')
% figure
% bar([spikeMeans{8,4},spikeMeans{10,4},spikeMeans{12,4}])
% title('OFF Smooth')
% 
% figure
% bar([spikeMeans{2,5},spikeMeans{4,5},spikeMeans{6,5}])
% figure
% bar([spikeMeans{8,5},spikeMeans{10,5},spikeMeans{12,5}])
% 
% figure
% bar([spikeMeans{2,6},spikeMeans{4,6},spikeMeans{6,6}])
% figure
% bar([spikeMeans{8,6},spikeMeans{10,6},spikeMeans{12,6}])
%% Kmeans sort

celltoAnalyze = ["OFFParasol","OFFSmooth"];

cellnameList2 = cellnameList;

kmeansGroups2 = kmeansGroups;
kmeansGroups2(contains(cellnameList2,celltoAnalyze)==0,:)=[];
cellnameList2(contains(cellnameList2,celltoAnalyze)==0)=[];
cellnameList2(sum(kmeansGroups2(:,1:3),2)<5)=[];
kmeansGroups2(sum(kmeansGroups2(:,1:3),2)<5,:)=[]
% kmeansGroups2(sum(kmeansGroups2(strcmp(cellnameList,'ONParasol')==0,1:3),2)>5,:)=[];
%% data edit after sorting into types 
parseNewData = 0;
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'Motion And Noise';
typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/');
recordingType='extracellular';

cd(typesFolder)
dataDir = dir;
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[]; %removes macOS hidden stuff .. sometimes 2 sometimes 3???
% folderSet(:,:,1:3)=[]; 
size(folderSet,3)

% splitByType{1,1}=["Strong ON OS"]
% splitByType{2,1}=["noiseClass","epochGroupLabel","frameDwell","apertureRadius","barOrientation","backgroundClass"]
splitFactors =["noiseClass","epochGroupLabel","frameDwell", "apertureRadius","barOrientation","backgroundClass"];
for typeIndex = 1:size(folderSet,3)
    
    cd(strcat(typesFolder,folderSet(:,:,typeIndex)))
    cellType = strtrim(folderSet(:,:,typeIndex));
    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    for cellIndex = 1:size(fileIndex,3)
    cellYear = fileIndex(:,1:4,cellIndex);
    cellMD = fileIndex(:,5:8,cellIndex);
    cellNum = fileIndex(:,9:11,cellIndex);
     
    experimentFolder = strcat([cellYear '_' cellMD '/']);
   

    load(strcat([exportFolder experimentFolder cellYear cellMD cellNum]))
    
    baseEpochs = epochs;
    
    load(strcat([exportFolder experimentFolder cellYear cellMD cellNum '_FT']))    
    
    frameEpochs = epochs;
    
    [splitCell,indexHolder,~,~,~] = makeData(baseEpochs,frameEpochs,protocolToAnalyze,splitFactors,6,recordingType,parseNewData);
%     [splitCell,indexHolder,spikingData,frameTimings,metaData] = makeData(spikeEpochs,frameTs,protocolToAnalyze,splitFactors,6,recordingType)
     save(strcat([typesFolder strtrim(folderSet(:,:,typeIndex)) '/' fileIndex(:,1:11,cellIndex),'_','MN','.mat']),'splitCell','indexHolder','-append')
    end
    
end
 


%% mistake eraser

%     cd(strcat(exportFolder,fileIndex(:,1:4,d),'_',fileIndex(:,5:8,d)))
    splitFactors = ["noiseClass","epochGroupLabel","frameDwell","backgroundClass"];
    
    load(strcat(fileIndex(:,1:11,d),'.mat'))
    spikeEpochs=epochs;
    load(strcat(fileIndex(:,1:11,d),'_FT.mat'))
    frameTs = epochs;
    saveLabel = 'MN';
    [splitCell,indexHolder,spikingData,frameTimings,metaData] = makeData(spikeEpochs,frameTs,protocolToAnalyze,splitFactors,6,recordingType);
    
%     cd(typeDir(1).folder)
    save(strcat(fileIndex(:,1:11,d),'_','MN','.mat'),'splitCell','indexHolder','frameTimings','-append')