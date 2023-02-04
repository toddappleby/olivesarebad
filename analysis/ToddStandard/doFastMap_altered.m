function [sta,dataAll,frameValues,frameValuesContinuous,frameTimesContinuous,staPerCheck,staPerCheck2] = doFastMap(r,numXChecks,numYChecks,numYStixels,numXStixels,numFrames,seeds,timings)

binRate = 1000;
binsPerFrame=1;
% Load the data.
% load('/Users/toddappleby/Downloads/ParasolOn_extracellular_20190115Ac6.mat');
% load('/Users/toddappleby/Downloads/SmoothMonoOff_extracellular_20190313Ac1.mat');

prePts = timings(1)*1e-3*binRate;
preFrames = timings(1)*1e-3*60;
numFramesOld = numFrames;
numFrames = numFrames*2;  %!!!
% numFrames2

frameTimes = round(cumsum([1 binRate/60*ones(1,preFrames*2+numFrames)]));
% frameTimes = round(cumsum([1 binRate/60*ones(1,preFrames*2+numFrames2)]));


frameTimes(frameTimes > (size(r,2)*2)) = []; %!!!

frameValuesContinuous=[];
frameTimesContinuous=[];
dataAll=[];
addFTs=[];

ft = zeros(numYChecks, numXChecks, numFramesOld*binsPerFrame);
size(frameTimes)
for m = 1 : length(seeds)
   
    
    [fv1,frameValues,~,fv2] = getFastNoiseFrames(numXStixels,numYStixels, numXChecks, numYChecks, 'BY', numFramesOld, 2, seeds(m));

    size(fv1);
    size(fv2);
    size(frameValues);
    
    
    if binsPerFrame > 1
        frameValues = upsampleFrames(frameValues,binsPerFrame);
    end
    % First frame.

    fstart = find(frameTimes > prePts,1);
    if fstart < 1, fstart = 1; end
    % Data binning per frame

    max(numFrames,length(frameTimes(fstart:end))-1)
    
    bdata = zeros(1,numFramesOld*binsPerFrame);
    
  
    for n = 2 : 2 : max(numFrames,length(frameTimes(fstart:end))-1)  %!!!
        idx = frameTimes(fstart+(n-1)) : frameTimes(fstart+n)-1;  
       
        for p = 1 : binsPerFrame  %!!!
            idx2 = floor((p-1)*length(idx)/binsPerFrame)+(1 : floor(length(idx)/binsPerFrame));
            bdata((n-1)*binsPerFrame+p) = mean(r(m,idx(idx2))); %!!!
        end
        
    end
    % revcor
    
    bdata = bdata - median(bdata(60*binsPerFrame+1:end));
    bdata(1:60*binsPerFrame) = 0;
    size(frameValues)
    size(bdata)

    for n = 1 : size(frameValues,1)
        for p = 1 : size(frameValues,2)
%             thing1=squeeze(ft(n,p,:)) + fft(bdata(:));
%             size(thing1)
%             size(conj(fft(squeeze(frameValues(:,n,p)))))

            ft(n,p,:) = squeeze(ft(n,p,:)) + fft(bdata(:)) ...
                .* conj(fft(squeeze(frameValues(n,p,:))));
            
            if n == 24 && p == 24
                tester = ft(n,p,:);
                staPerCheck(m,:) =ifft(tester);
            elseif n == 24 && p == 23
                tester = ft(n,p,:);
                staPerCheck2(m,:) =ifft(tester);
            end
            
        end 
    end
  

%     tester=ft(24,24,:);
%     tester=ifft(tester(:));
%     plot(tester)
%   pause
    
%     staEachSeed = zeros(numYChecks,numXChecks,30);
%     for k = 1 : numYChecks
%         for m = 1 : numXChecks
%             tmp = ifft(squeeze(ft(k,m,:)));
% %              plot(tmp(1:30))  
% %              pause
%             staEachSeed(k,m,:) = tmp(1 : 30);
%         end
%     end
%     
% m

% frameTimes(frameTimes<250)=[];
% frameTimes(frameTimes>21250) = [];


%   frameValuesContinuous(:,:,(1+numFrames*(m-1)):(numFrames*m)) = frameValues;

 
%   dataAll = [dataAll bdata]; %#ok<AGROW>

% size((frameTimes+(frameTimes*(m-1))))
%   if isempty(frameTimesContinuous)
%   frameTimesContinuous = [frameTimesContinuous (frameTimes+(frameTimes*(m-1)))]; %#ok<AGROW>
%   else
%       size(frameTimesContinuous)
%       m
%   frameTimesContinuous = [frameTimesContinuous (frameTimesContinuous+(frameTimes+(frameTimes*(m-1))))]; %#ok<AGROW>
%   end
% frameTimes+((size(frameTimes,2)*(m-1)))
% addFTs = frameTimes + (max(frameTimes)*(m-1));

% frameTimesContinuous=cat(2,frameTimesContinuous,addFTs);


end

sta = zeros(numYChecks,numXChecks,30);
    for k = 1 : numYChecks
        for m = 1 : numXChecks
            tmp = ifft(squeeze(ft(k,m,:)));
%              plot(tmp(1:30))  
%              pause
            sta(k,m,:) = tmp(1 : 30);
        end
    end

% image code
imagesc(sta(:,:,1))
colormap(parula)
colorbar
title('Linear Prediction at each Spatial Noise Frame','fontSize',17)
xlabel('X Frame Coordinates','fontSize',14)
ylabel('Y Frame Coordinates','fontSize',14)

for g = 1:4
subplot(2,2,g)
imagesc(sta(:,:,g+2))
pause
colormap(parula)
colorbar
end

frameTimesContinuous = []; % comment out if doing something in python or w/ a CNN again

end
