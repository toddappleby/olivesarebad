function centeringBars(XspikeMatrix,YspikeMatrix,positionX,positionY,timings,tFrequency,sampleRate)

% sampleRate = 10000;
% 
% stimStart = (timings(1)*1e-3)*sampleRate+1;
% stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
% stimOff = (timings(1)+timings(2)+timings(3))*10;
% 
%             data = mean(XspikeMatrix(sortedIndex(finalInd),:),1);
%                     binnedData = BinSpikeRate(data(stimStart:stimEnd), 100, sampleRate);
%                 [F, phase] = frequencyModulation(binnedData, ...
%                 100, indexHolder{1,a}(1,params.tfreq(1)), 'avg', 1:2, []);
%                 avgF1(u)= F(1);
%                 avgF2(u)= F(2);    
                


sender = struct();
uniqueP = unique(positionX);
uniqueP = sort(uniqueP);
for s = 1:length(uniqueP)
barIndex = find(positionX == uniqueP(s));
barIndex = barIndex';
sender(s).data = mean(XspikeMatrix(barIndex,:),1);
sender(s).params.spatialFrequency = uniqueP(s);
sender(s).params.temporalFrequency = tFrequency;
sender(s).params.preTime = timings(1)*1e-3;
sender(s).params.stimStart = (timings(1)*1e-3)*sampleRate+1;
sender(s).params.stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
sender(s).params.sampleRate = sampleRate;
end
result = sMTFAnalysisTodd(sender,1e4,'spikes','avg');
figure(1)
subplot(1,2,1)
    plot(result.uniqueSF,result.avgF1);
    hold on
    plot(result.uniqueSF,result.avgF2);
    title('F1 Response to Spot')
    xlabel('Spot Radius (micron)')
    ylabel('F1 Amplitude (Hz)')

uniqueP = unique(positionY);
uniqueP = sort(uniqueP);
    
for s = 1:length(uniqueP)
    barIndex = find(positionY == uniqueP(s));
    barIndex = barIndex';
    sender(s).data = mean(YspikeMatrix(barIndex,:),1);

    sender(s).params.spatialFrequency = uniqueP(s);
    sender(s).params.temporalFrequency = tFrequency;
    sender(s).params.preTime = timings(1)*1e-3;
    sender(s).params.stimStart = (timings(1)*1e-3)*sampleRate+1;
    sender(s).params.stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
    sender(s).params.sampleRate = sampleRate;
end

result = sMTFAnalysisTodd(sender,1e4,'spikes','avg');
subplot(1,2,2)
    plot(result.uniqueSF,result.avgF1);
    hold on
    plot(result.uniqueSF,result.avgF2);
    title('F1 Response to Spot')
    xlabel('Spot Radius (micron)')
    ylabel('F1 Amplitude (Hz)')

end

