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
figure(3)
subplot(13,1,s)
plot(psth(XspikeMatrix(barIndex,:),6+2/3,10000,1))
sgtitle('x bar psth at each bar position')

ylabel(uniqueP(s))


end
result = sMTFAnalysisTodd(sender,1e4,'spikes','avg');
figure(1)
subplot(2,2,1)
    plot(result.uniqueSF,result.avgF1);
    hold on
    plot(result.uniqueSF,result.avgF2);
    title('Response to 1D Bar (X Position)')
    xlabel('Vertical Bar Position (micron)')
    ylabel('F1 (Blue) & F2 (Red) Amplitude (Hz)')
    
subplot(2,2,3)
plot(result.uniqueSF,result.avgPh)
  xlabel('Vertical Bar Position (micron)')
    ylabel('Phase Amplitude')


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
    figure(4)
    subplot(13,1,s)
    plot(psth(YspikeMatrix(barIndex,:),6+2/3,10000,1))
    sgtitle('y bar psth at each bar position')
    ylabel(uniqueP(s))
end

result = sMTFAnalysisTodd(sender,1e4,'spikes','avg');
figure(1)
subplot(2,2,2)
    plot(result.uniqueSF,result.avgF1);
    hold on
    plot(result.uniqueSF,result.avgF2);
    title('Response to 1D Bar (Y Position)')
    xlabel('Horizontal Bar Position (micron)')
    ylabel('F1 (Blue) & F2 (Red) Amplitude (Hz)')
    
subplot(2,2,4)
plot(result.uniqueSF,result.avgPh)
    xlabel('Horizontal Bar Position (micron)')
    ylabel('Phase Amplitude')
end

