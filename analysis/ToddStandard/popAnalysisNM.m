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
recordingType = 'extracellular';
binRate = 1e3;
% size(folderSet,3)
for c = 2:2
    
    cd(strcat(typesFolder,folderSet(:,:,c)))
    
    cellType = strtrim(folderSet(:,:,c))
    
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == '2');
    
    cellTypeData = [];
    

        for d = 1:size(fileIndex,3)
           outI = []; 

         load(strtrim(fileIndex(:,:,d)))
         %Where params are most appropriate?  prob the ones I ran the most:
         clear indLength
            for e = 1: length(indexHolder)
                indLength(e) = length(indexHolder{2,e});
            end
            
            maxI = maxk(indLength,3); %top 3 conditions -- must be best options for seq rand & static
            logI = ismember(string(indLength),string(maxI));
            outI = find(logI==1); %index holder indices of top options
            
            %which indices correspond to which condition?
            clear holdLabel
            for ee = 1:length(outI)
                holdLabel(ee) = string(indexHolder{1,outI(ee)}(4));
            end
            
            seqIndex = indexHolder{2,outI(strcmp(holdLabel,"sequential"))};
            randIndex = indexHolder{2,outI(strcmp(holdLabel,"random"))};
            staticIndex = indexHolder{2,outI(strcmp(holdLabel,"stationary"))};
            
            
            
            prePts = 250 * 1e-3 * binRate;
            stimPts = 10000 * 1e-3 * binRate;
            tailPts = 250 * 1e-3 * binRate;
           lfilter=[];
            for f = 1:3
                indices = indexHolder{2,outI(f)};
                frameDwell = ones(length(indices),1)*double(indexHolder{1,outI(f)}(3));
                
                noiseVars = struct();
                noiseVars.contrast = 0.3333;
                noiseVars.type = indexHolder{1,outI(f)}(1);
                
                
                frames = manookinlab.ovation.getFrameTimesFromMonitor(frameTimings(indices+size(spikingData,1),:), metaData.sampleRate, binRate);
                frameSeq = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),binRate,frames,1,seed(indices),frameDwell);
                
                response = zeros(size(indices,1),10500);
                for ff = 1:length(indices)
                    response(ff,:) = binSpikeCount(spikingData(indices(ff),:)/metaData.sampleRate, binRate, metaData.sampleRate);
                    response(ff,:) = psth(response(ff,:)*binRate,6+2/3,binRate,1);
                end
                
                
             stimulus = frameSeq(seed(indices) ~= 1, :);
             responses = response(seed(indices) ~= 1, :);
                
            lfilter(f,:) = getLinearFilter(stimulus, responses, ...
                'analysisType', 'revcorr', ...
                'fourierCorrection',false, ...
                'binRate', binRate, ...
                'filterTime', 0.5, ...
                'frameRate', metaData.frameRate);
            lfilter(f,:) = lfilter(f,:) / norm(lfilter(f,:));
            
            stimIndex = prePts + (1 : stimPts);
            
            
                
            end
            
            lfilter=mean(lfilter,1);
   
            
            %doing it twice because no time to rewrite simple thing?
            for g = 1:3
                
                                indices = indexHolder{2,outI(g)};
                frameDwell = ones(length(indices),1)*double(indexHolder{1,outI(g)}(3));
                
                noiseVars = struct();
                noiseVars.contrast = 0.3333;
                noiseVars.type = indexHolder{1,outI(g)}(1);
                
                
                frames = manookinlab.ovation.getFrameTimesFromMonitor(frameTimings(indices+size(spikingData,1),:), metaData.sampleRate, binRate);
                frameSeq = getTemporalNoiseFramesClarinet(noiseVars,timings(1),timings(2),timings(3),binRate,frames,1,seed(indices),frameDwell);
                
                response = zeros(size(indices,1),10500);
                for gg = 1:length(indices)
                    response(gg,:) = binSpikeCount(spikingData(indices(gg),:), binRate, metaData.sampleRate);
                    response(gg,:) = psth(response(gg,:),6+2/3,binRate,1);
                end
                
                S = frameSeq;
                R = response;
                lf = lfilter(1 : size(S,2));
                P = zeros(size(R));
                for p = 1 : size(S,1)
                    tmp = ifft( fft(S(p,:)) .* fft(lf) );
                    P(p,:) = real(tmp);
                    P(p,:) = P(p,:)./std(P(p,:));
                end
                [xBin(g,:),yBin(g,:)] = binNonlinearity(P(:,stimIndex(501:end)),R(:,stimIndex(501:end)),nonlinearityBins);
                
                
%                 nlParams(g,:) = models.ln.fitNonlinearityParams(xBin(g,:), yBin(g,:));
%                 prediction(indices,:) = P;
            end
                
            
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
              
              
              MSEdecrease(d) = (MSEseqGain-MSEseqHZ)/MSEseqGain;
              
              
              
              %both 
%                 % rand seq
%               mParams = fitMultiVarParams(xBin([1,2],:),yBin([1,2],:),2,starterParams);
%               outNLHZGain = multiVarNL(mParams,xBin([1,2],:));

                % rand static
              mParams = fitMultiVarParams(xBin([1,3],:)',yBin([1,3],:)',2,starterParams);
              outNLBoth = multiVarNL(mParams,xBin([1,3],:)');
              
            RANDpChangeHZ(d)  = (mParams(5)-mParams(3))/mParams(5);
              RANDpChangeGain(d) = (mParams(4)-mParams(2))/(-1*mParams(4));
              
              
                % seq static
              mParams = fitMultiVarParams(xBin([2,3],:)',yBin([2,3],:)',2,starterParams);
              outNLBoth = multiVarNL(mParams,xBin([2,3],:)');
              
            
             SEQpChangeHZ(d)  = (mParams(5)-mParams(3))/mParams(5);
              SEQpChangeGain(d) = (mParams(4)-mParams(2))/(-1*mParams(4));
              %percent change final-initial/initial
            
        end
        
        RandHZOut=RANDpChangeHZ(1:d)'*100;
        RandGainOut = RANDpChangeGain(1:d)'*100;
        SeqHZOut = SEQpChangeHZ(1:d)'*100;
        SeqGainOut = SEQpChangeGain(1:d)'*100;
        meanDecrease = mean(MSEdecrease)*100;
        errorDecrease = sem(MSEdecrease)*100;
        
        
        save(strcat(cellType,'_MN'),'RandHZOut','RandGainOut','SeqHZOut','SeqGainOut','meanDecrease','errorDecrease')
        
end     

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
    
%     cd(typeDir(1).folder)
    save(strcat(fileIndex(:,1:11,d),'_','MN','.mat'),'splitCell','indexHolder','frameTimings','-append')
