function simpleResponse(spikeMatrix,psthMatrix,epochStorage,timings,splitCell,indexHolder,params)

sampleRate = 10000;
stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;
a=0;
numSubPlots = 3;

% figure(11)
% subplot(2,1,1)

xvals = linspace(0,stimOff/10,length(psthMatrix));


for b = 1:size(indexHolder,2)
    
    
    xVar = splitCell{2,size(splitCell,2)};
    xVarStatic = unique(xVar);
    xVar = xVar(indexHolder{2,b});
    [xVar, I] = sort(xVar);
    sortedIndex = indexHolder{2,b};
    sortedIndex = sortedIndex(I);
     
    if size(splitCell,2) == 1
        xVar = splitCell{2,size(splitCell,2)};
    end
    
    spikevalsforFit = []; 
    

    for g = 1:length(xVarStatic)
        
        if isnan(unique(str2double(xVar)))
            xVarLabel = unique(xVar);
        else
            xVarLabel = unique(str2double(xVar));
        end

        finalInd = find(xVar==string(xVarLabel(g)));
        
        if size(splitCell,2) == 1 %makes this work w/ single splitter.  indices in index holder don't need separating by sorting var (like orientation)
            sortedIndex = indexHolder{2,g};
            finalInd = 1:1:length(sortedIndex);
        end
            
        
        
        
        
    figure(1+b)
    
    subplot(length(xVarLabel),1,g)
    plot(xvals,mean(psthMatrix(sortedIndex(finalInd),:),1));
    xline(stimStart/10,'LineStyle','--','Color','r','LineWidth',2);
    xline(stimEnd/10,'LineStyle','--','Color','r','LineWidth',2);
    ylabel(xVarLabel(g));
    xlabel('time (ms)');
    sgtitle(strcat([splitCell{1,1:size(splitCell,2)-1}]," ",indexHolder{1,b}));
    
    figure(20+g+(b*3))
    plot(xvals,mean(psthMatrix(sortedIndex(finalInd),:),1));
    

    figure(15+a)  

    set(gcf,'Position', [100, 100, 1024, 1200])
    % spikedThing = sum(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1)/max(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1)),2);
    subplot(numSubPlots,1,b-(numSubPlots*a))
    spikedThing = sum(mean(spikeMatrix(sortedIndex(finalInd),(params.moveTime+500):end),1));
    spikeError = sem(sum(spikeMatrix(sortedIndex(finalInd),(params.moveTime+500):end),2));
    scatter(g,spikedThing,50,'filled')
    errorbar(g,spikedThing,spikeError,spikeError)
    xticks(1:1:length(xVarLabel))
    xticklabels(xVarLabel);
    ylabel('mean spike count')
    title(strcat([splitCell{1,1:size(splitCell,2)-1}]," ",indexHolder{1,b}));
    sgtitle(strcat("Spike Count for each: ",splitCell{1,size(splitCell,2)}))
    
    spikevalsforFit = [spikevalsforFit spikedThing];
    
    hold on
    
spikedThing2(g,:) = spikedThing;
spikedError2(g,:) = spikeError;
    
    end
    
    outSpike = spikedThing2';
    outSpike = outSpike/max(outSpike);
    outError = spikedError2';
    
    save('OMDfile','outSpike','outError')
    
   figure(15+a)
   pVals = polyfit(xticks,spikevalsforFit,4);
   pFit = polyval(pVals,xticks);
  
   plot(xticks,pFit,'LineStyle','--','LineWidth',.5)
   
    if mod(b,3) == 0
        a = a +1;
        
    end
   
    
end
figure
if params.saveGraph == 1 && strcmp(params.protocolID,"Single Spot") || strcmp(params.protocolID,"Chromatic Spot")
    for L = 1:size(epochStorage,1)
    
    plot(epochStorage(L,:)+(L*300))
    hold on

    end
    epochStorage = epochStorage';
    save('epochStorage.mat','epochStorage')
end
    

%
% if strcmp(params.protocolID,'Single Spot')
% 
% spotThing = mean(psthMatrix(indexHolder{2,1},:),1);
% spotThing = spotThing';
% % saveName = strcat('spotResponse_',params.cellName);
% save('spotResponse.mat','spotThing')
% elseif strcmp(params.protocolID,'Led Pulse')
%     
% uvResponse = mean(psthMatrix(indexHolder{2,1},:),1);
% uvResponse = uvResponse';
% 
% redResponse = mean(psthMatrix(indexHolder{2,2},:),1);
% redResponse = redResponse';   
% save('ledResponse.mat','uvResponse','redResponse')
% else
%     
% uvResponse = mean(psthMatrix(indexHolder{2,1},:),1);
% uvResponse = uvResponse';
% 
% redResponse = mean(psthMatrix(indexHolder{2,2},:),1);
% redResponse = redResponse';    
% save('chromaticResponse.mat','uvResponse','redResponse')
% end
%ORIGINAL SCATTER:
%  
% for s = 1:size(splitCell,2)
%     for t = 1:length(indexHolder)
% 
% 
%     figure(10+s)   
%     % spikedThing = sum(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1)/max(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1)),2);
%     spikedThing = sum(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1));
%     scatter(t,spikedThing,50,'filled')
%     xticks([1:1:length(indexHolder)])
%     labelHolder(1,t)=indexHolder{1,t};
%     hold on
% 
%     end
%     xticklabels(labelHolder);
% end



% plot(xvals,mean(psthMatrix(indexHolder{2,2},:)))
% xline(stimStart/10,'LineStyle','--','Color','m','LineWidth',2)
% xline(stimEnd/10,'LineStyle','--','Color','m','LineWidth',2)
% xlabel('time (ms)')
% ylabel(indexHolder{1,2})
% subplot(2,1,2)
% plot(xvals,mean(psthMatrix(indexHolder{2,3},:)))
% xline(stimStart/10,'LineStyle','--','Color','m','LineWidth',2)
% xline(stimEnd/10,'LineStyle','--','Color','m','LineWidth',2)
% xlabel('time (ms)')
% ylabel(indexHolder{1,3})


end