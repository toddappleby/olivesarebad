clear indLength 
binRate=1e3;
timings=[250,10000,250];
nonlinearityBins=100;
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
            holdLabel = holdLabel(1:3);
            seqIndex = indexHolder{2,outI(strcmp(holdLabel,"sequential"))};
            randIndex = indexHolder{2,outI(strcmp(holdLabel,"random"))};
            staticIndex = indexHolder{2,outI(strcmp(holdLabel,"stationary"))};
            
            prePts = 250 * 1e-3 * binRate;
            stimPts = 10000 * 1e-3 * binRate;
            tailPts = 250 * 1e-3 * binRate;
           lfilter=[]
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
                    P(p,:) = P(p,:)/std(P(p,:));
                end
                [xBin(g,:),yBin(g,:)] = binNonlinearity(P(:,stimIndex(501:end)),R(:,stimIndex(501:end)),nonlinearityBins);
                
                
                nlParams(g,:) = models.ln.fitNonlinearityParams(xBin(g,:), yBin(g,:));
                prediction(indices,:) = P;
            end
            
figure(10)
plot(xBin(strcmp(holdLabel,'random'),:),yBin(strcmp(holdLabel,'random'),:),'b')
hold on

plot(xBin(strcmp(holdLabel,'sequential'),:),yBin(strcmp(holdLabel,'sequential'),:),'r')
plot(xBin(strcmp(holdLabel,'stationary'),:),yBin(strcmp(holdLabel,'stationary'),:),'k')

xNorm(1,:) = xBin(strcmp(holdLabel,'random'),:);
xNorm(2,:) = xBin(strcmp(holdLabel,'sequential'),:);
xNorm(3,:) = xBin(strcmp(holdLabel,'stationary'),:);

xNorm = xNorm';

WHY(1,:)=yBin(strcmp(holdLabel,'random'),:);
WHY(2,:)=yBin(strcmp(holdLabel,'sequential'),:);
WHY(3,:)=yBin(strcmp(holdLabel,'stationary'),:);

WHY=WHY';


starterParams = [1,-1.5];
mParams = fitMultiVarParams(xBin([1,2],:)',yBin([1,2],:)',1,starterParams);
outNLHZ = multiHZNL(mParams,xBin([1,2],:)');

 MSErandHZ = immse(yBin(1,:)',outNLHZ(:,1))
  MSEseqHZ = immse(yBin(2,:)',outNLHZ(:,2))

figure(10)
hold on
plot(xBin(1,:),outNLHZ(:,1),'--')
plot(xBin(2,:),outNLHZ(:,2),'--')

xBin1 = xBin';
outNLHZ1 = outNLHZ';

tester = mParams;

figure(11)
plot(xBin(strcmp(holdLabel,'random'),:),yBin(strcmp(holdLabel,'random'),:),'b')
hold on

plot(xBin(strcmp(holdLabel,'sequential'),:),yBin(strcmp(holdLabel,'sequential'),:),'r')
plot(xBin(strcmp(holdLabel,'stationary'),:),yBin(strcmp(holdLabel,'stationary'),:),'k')

starterParams = [.2,-1.5];
mParams = fitMultiVarParams(xBin([1,2],:)',yBin([1,2],:)',0,starterParams);
outNLGain = multiGainNL(mParams,xBin([1,2],:)');

 MSErandGain = immse(yBin(1,:)',outNLGain(:,1))
  MSEseqGain = immse(yBin(2,:)',outNLGain(:,2))

figure(11)
hold on
plot(xBin(1,:),outNLGain(:,1),'--')
plot(xBin(2,:),outNLGain(:,2),'--')

outNLGain2 = outNLGain';
xBin2 = xBin';

save('bigSaver','xNorm','WHY','xBin1','outNLHZ','xBin2','outNLGain');

%%

fredDots(1,:)=[1 .83 .58 .54 .35 .42 .4 .43];

fredDots(2,:)=[1 .88 .82 .69 .65 .59 .61 .69]

fredDots(3,:)=[1 .92 .75 .53 .33 .34 .3 .23]

fredDots(4,:)=[1
0.7999999999999999
0.6692307692307692
0.4769230769230769
0.5307692307692308
0.5076923076923077
0.4153846153846154
0.5461538461538461]


