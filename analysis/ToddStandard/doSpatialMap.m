function [sta,frameValues] = doSpatialMap(r,numXChecks,numYChecks,numFrames,seeds,timings)

binRate = 1000;
binsPerFrame=1;
% Load the data.
% load('/Users/toddappleby/Downloads/ParasolOn_extracellular_20190115Ac6.mat');
% load('/Users/toddappleby/Downloads/SmoothMonoOff_extracellular_20190313Ac1.mat');

prePts = timings(1)*1e-3*binRate;
preFrames = timings(1)*1e-3*60;
frameTimes = round(cumsum([1 binRate/60*ones(1,preFrames*2+numFrames)]));
frameTimes(frameTimes > size(r,2)) = [];

% ft = zeros(numYChecks, numXChecks, numFrames*binsPerFrame+50);
ft = zeros(numYChecks, numXChecks, numFrames*binsPerFrame);
for m = 1 : length(seeds)
    
    frameValues = getSpatialNoiseFrames(numXChecks, numYChecks, numFrames, 'binary', 'achromatic', seeds(m));
%     frameValues = getJitteredNoiseFrames(24, 24, 46, 46, numFrames, 2, seeds(m));
%     [fv1,frameValues,~,fv2] = getFastNoiseFrames(26,17, 50, 32, 'BY', numFrames, 2, seeds(m));
    
    if binsPerFrame > 1
        frameValues = upsampleFrames(frameValues,binsPerFrame);
    end
    % Find the start frame.
    fstart = find(frameTimes > prePts,1);
    if fstart < 1, fstart = 1; end
    % Bin the data.
    bdata = zeros(1,numFrames*binsPerFrame);
    for n = 1 : min(numFrames,length(frameTimes(fstart:end))-1)
        idx = frameTimes(fstart+(n-1)) : frameTimes(fstart+n)-1;
        for p = 1 : binsPerFrame
            idx2 = floor((p-1)*length(idx)/binsPerFrame)+(1 : floor(length(idx)/binsPerFrame));
            bdata((n-1)*binsPerFrame+p) = mean(r(m,idx(idx2)));
        end
    end
    % Do reverse correlation.
    bdata = bdata - median(bdata(60*binsPerFrame+1:end));
    bdata(1:60*binsPerFrame) = 0;

    for n = 1 : size(frameValues,2)
        for p = 1 : size(frameValues,3)
               ft(n,p,:) = squeeze(ft(n,p,:)) + fft([zeros(25,1);bdata(:);zeros(25,1)]) ...
                .* conj(fft([zeros(25,1);squeeze(frameValues(:,n,p));zeros(25,1)]));
        end
    end
end

sta = zeros(numYChecks,numXChecks,30);
    for k = 1 : numYChecks
        for m = 1 : numXChecks
            tmp = ifft(squeeze(ft(k,m,:)));
            sta(k,m,:) = tmp(1 : 30);
        end
    end
% lobePts = round(0.05*30/0.5) : round(0.15*30/0.5);
% spatialRFx = squeeze(mean(filterImage(:,:,lobePts),3));
% image code
imagesc(sta(:,:,1))
colormap(parula)
colorbar
title('Linear Prediction at each Spatial Noise Frame','fontSize',17)
xlabel('X Frame Coordinates','fontSize',14)
ylabel('Y Frame Coordinates','fontSize',14)
end
