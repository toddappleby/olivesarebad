function frameSeq = getTemporalNoiseFramesClarinet(noiseVars, preTime, stimTime, tailTime, binRate, frameTimes, frameShift,seedList,frameDwell)


if ~exist('frameShift','var')
    frameShift = 0;
end

prePts = round(preTime * 1e-3 * binRate);
stimPts = round(stimTime * 1e-3 * binRate);
samplePts = floor((preTime + stimTime + tailTime) * 1e-3 * binRate);


frameSeq = zeros(length(seedList), samplePts);
for k = 1 : length(seedList)
    
    noiseContrast = noiseVars.contrast;
    noiseClass = noiseVars.type;
    seed = seedList(k);
    
    % Seed the random number generator.
    noiseStream = RandStream('mt19937ar', 'Seed', seed);
    
    % Reconstruct the frame sequence.
    fTimes = frameTimes{k};
   
    firstFrame = find(fTimes > prePts,1) - frameShift;
    lastFrame = find(fTimes > prePts+stimPts,1);
    
    fTimes = fTimes(firstFrame : lastFrame);
    numFrames = length(fTimes);
    
    nframes = length(fTimes);
    
    switch noiseClass
        case 'gaussian'
            fSeq = 0.3*noiseStream.randn(1,nframes);
        case 'binary'
            fSeq = 2*(noiseStream.rand(1,nframes)>0.5)-1;
        case 'uniform'
            fSeq = 2*(noiseStream.rand(1,nframes))-1;
    end
    
   
    fSeq = noiseContrast * fSeq;
    
    if frameDwell(k) > 1
     fSeq = ones(frameDwell(k),1)*fSeq(:)';
    end
    
    fSeq(fSeq > 1) = 1;
    fSeq(fSeq < -1) = -1;
    
    for m = 1 : numFrames-1
        if binRate > 60
            frameSeq(k, fTimes(m) : fTimes(m+1)) = fSeq(m);
        else
            frameSeq(k, m+prePts-frameShift) = fSeq(m);
        end
    end
end
