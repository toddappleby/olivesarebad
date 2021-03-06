function simpleSpots(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

sampleRate = 10000;
stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;

% to indicate pos or neg contrast spot for each protocol:
if strcmp(params.stimName,"mTF")
analysisSplitters = string(splitCell(1,:));
spotIntensityPlace = find(strcmp(analysisSplitters,"contrast")==1);
else
analysisSplitters = string(splitCell(1,:));
spotIntensityPlace = find(strcmp(analysisSplitters,"spotIntensity")==1);
end

    for a = 1:size(indexHolder,2)
        spotSizes = splitCell{2,size(splitCell,2)};
        spotSizes = spotSizes((indexHolder{2,a}));
         [spotSizes I] = sort(spotSizes);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);  %sort indices like spotSizes (which are low to high)
         
         
         for b = 1:length(unique(spotSizes))
            uniqueSpotSizes = unique(spotSizes); %because I don't know how to index the unique function
            finalInd = find(spotSizes==uniqueSpotSizes(b));    
            onResp(b) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
            offResp(b) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));                
            
            if params.dataType ==1
            psthData(b,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(b,:)+(100*(b-1)))
            title(indexHolder{1,a})
            else
            psthData(b,:) = mean(spikeMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(b,:)+(1*(b-1)))
            title(indexHolder{1,a})
            end
            
            
            
         end
        
         
         
         figure(1) 
        
         subplot(size(indexHolder,2),1,a)
         if indexHolder{1,a}(spotIntensityPlace) <= 0
        plot(uniqueSpotSizes,onResp,'Color','r','LineWidth',2); hold on
        plot(uniqueSpotSizes,offResp,'Color','b','LineWidth',2);
        legend('Off Response','On Response')
        xlabel('spot sizes')
        ylabel('spike count')
        title('Dark Spot')
%         saveName = strcat('expandingDarkSpot_',params.cellName);
        onRespFlip = offResp';
        offRespFlip = onResp';
        uniqueSpotSizesFlip = uniqueSpotSizes';
        save('expandingDarkSpots.mat','uniqueSpotSizesFlip','onRespFlip','offRespFlip')
        
        %dog
%        size(onResp)
%        size(unique(spotSizes))
%         [Kc, sigmaC, Ks, sigmaS, baseFiring] = fitDoGAreaSummation2(unique(spotSizes),onResp,0);
%         
%         fitX = 0:max(spotSizes);
%         fitY = DoGAreaSummation([Kc sigmaC Ks sigmaS baseFiring],fitX);
%         
%         figure(99)
%         plot(fitX,fitY);
        
         else
        plot(uniqueSpotSizes,onResp,'Color','b','LineWidth',2); hold on
        plot(uniqueSpotSizes,offResp,'Color','r','LineWidth',2);
        legend('On Response','Off Response')
        xlabel('spot sizes')
        ylabel('spike count')
        
      
        
        title("Light Spot --" + indexHolder{1,a}(2)) %will need to automate later
       
%         saveName = strcat('expandingLightSpot_',params.cellName);
        onRespFlip = onResp';
        offRespFlip = offResp';
        uniqueSpotSizesFlip = uniqueSpotSizes';
        save('expandingLightSpots.mat','uniqueSpotSizesFlip','onRespFlip','offRespFlip')
        
         end
      
    end
         

end