function freqAnalysis(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)
sampleRate = 10000;

stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;


analysisSplitters = string(splitCell(1,:));
tFreq= find(strcmp(analysisSplitters,"temporalFrequency")==1);
figure(10); clf; hold on


    for a = 1:size(indexHolder,2)
        radii = splitCell{2,size(splitCell,2)};
        radii = radii((indexHolder{2,a}));
         [radii I] = sort(radii);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
         
         

        for u = 1:length(unique(radii))
            uniqueRadii = unique(radii); %because I don't know how to index the unique function
            finalInd = find(radii==uniqueRadii(u));
%             angleIndex = find(orientation == uniqueAngle(u));
%             masterIndex = intersect(angleIndex,apertureIndex(a,:));
            data = mean(spikeMatrix(sortedIndex(finalInd),:),1);
                    binnedData = BinSpikeRate(data(stimStart:stimEnd), 100, sampleRate);
                [F, phase] = frequencyModulation(binnedData, ...
                100, indexHolder{1,a}(1,params.tfreq(1)), 'avg', 1:2, []);
                avgF1(u)= F(1);
                avgF2(u)= F(2) ;   
                
                phaser(u,:) = phase;
                spikeCount(u) = sum(data(stimStart:stimEnd),2);
                
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(16); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
        end
        
%         uniqueRadii

    figure(a)
    subplot(3,1,1);
    plot(uniqueRadii, avgF1(1:length(uniqueRadii)));
%     titleString = strcat({'F1 Resp'}, num2str(indexHolder{1,1}));
    title('F1 Resp')
    xlabel(analysisSplitters(1:3))
    hold on
    subplot(3,1,2);
    plot(uniqueRadii, avgF2(1:length(uniqueRadii)));
%     titleString = strcat({'F2 Resp'}, num2str(indexHolder{1,1}));
    title('F2 Resp')
    xlabel(num2str(indexHolder{1,a}'))
    
    subplot(3,1,3);
    plot(uniqueRadii, phaser(1:length(uniqueRadii)));
%     titleString = strcat({'Spike Count'}, num2str(indexHolder{1,1}));
    title('p h a s e')
    xlabel(analysisSplitters(1:3))
    
   figure(a+3)
   subplot(2,1,1)
   plot(uniqueRadii, avgF1(1:length(uniqueRadii)),'Color','b');
   hold on 
   plot(uniqueRadii, avgF2(1:length(uniqueRadii)),'Color','r');
   subplot(2,1,2)
   plot(uniqueRadii, phaser(1:length(uniqueRadii)));
   
%    saveName = strcat('mtfSpot_',params.cellName);
   savedF1 = avgF1(1:length(uniqueRadii));
   savedF1=savedF1';
   uniqueRadiiSave = uniqueRadii';
   save('mtfSpot.mat','uniqueRadiiSave','savedF1')
   
   
   
    
    end

end