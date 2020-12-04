function orientedStim(epochStorage,spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

% for z = 1:length(splitCell)-1    
%     uniqueSplit = unique(splitCell{2,z}); 
%     for c = 1:length(uniqueSplit)
%     splitIndex = find(splitCell{2,z}==uniqueSplit(c));
%     
%     end
% end
sampleRate = 10000;

stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;


analysisSplitters = string(splitCell(1,:));
tFreq= find(strcmp(analysisSplitters,"temporalFrequency")==1);
figure(10); clf; hold on
if strcmp(params.stimName, 'Grating')

    for a = 1:size(indexHolder,2)
        orientations = splitCell{2,size(splitCell,2)};
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
         
         

        for u = 1:length(unique(orientations))
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
%             angleIndex = find(orientation == uniqueAngle(u));
%             masterIndex = intersect(angleIndex,apertureIndex(a,:));
            data = mean(spikeMatrix(sortedIndex(finalInd),:),1);
                    binnedData = BinSpikeRate(data(stimStart:stimEnd), 100, sampleRate);
                [F, phase] = frequencyModulation(binnedData, ...
                100, str2num(indexHolder{1,a}(1,tFreq)), 'avg', 1:2, params.killCycle1);
                avgF1(u)= F(1);
                avgF2(u)= F(2);     
                spikeCount(u) = sum(data(stimStart:stimEnd),2);
                
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
           
            figure(15); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
        end

        
        uniqueOrientations = sort(str2double(uniqueOrientations));
        
    figure(a)
    subplot(1,3,1);
    polar(uniqueOrientations * 2 * pi/360, avgF1(1:length(uniqueOrientations))*1e-4);
%     titleString = strcat({'F1 Resp'}, num2str(indexHolder{1,1}));
    title('F1 Resp')
    xlabel(analysisSplitters(1:3))
    hold on
    subplot(1,3,2);
    polar(uniqueOrientations * 2 * pi/360, avgF2(1:length(uniqueOrientations))*1e-4);
%     titleString = strcat({'F2 Resp'}, num2str(indexHolder{1,1}));
    title('F2 Resp')
    xlabel(indexHolder{1,a}')
    
    subplot(1,3,3);
    polar(uniqueOrientations * 2 * pi/360, spikeCount(1:length(uniqueOrientations)));
%     titleString = strcat({'Spike Count'}, num2str(indexHolder{1,1}));
    title('Spike Count')
    xlabel(analysisSplitters(1:3))
    
    figure(10)

   figure(11)
   polar(uniqueOrientations * 2 * pi/360, avgF1(1:length(uniqueOrientations))*1e-4/max(avgF1(1:length(uniqueOrientations))*1e-4))
   hold on
   polar(uniqueOrientations * 2 * pi/360, avgF2(1:length(uniqueOrientations))*1e-4/max(avgF2(1:length(uniqueOrientations))*1e-4));
   polar(uniqueOrientations * 2 * pi/360, spikeCount(1:length(uniqueOrientations))/max(spikeCount(1:length(uniqueOrientations))));
    title('normalized F1,F2,Spikes')
    legend('F1','F2','Spikes','Location','NorthEastOutside')
    if a == 1
    normF1 = avgF1(1:length(uniqueOrientations))*1e-4/max(avgF1(1:length(uniqueOrientations))*1e-4);
    normF2 = avgF2(1:length(uniqueOrientations))*1e-4/max(avgF2(1:length(uniqueOrientations))*1e-4);
    normSpikes = spikeCount(1:length(uniqueOrientations))/max(spikeCount(1:length(uniqueOrientations)));
    normF1 = normF1';
    normF2 = normF2';
    normSpikes = normSpikes';
%     saveName = strcat('Gratings_',params.cellName);
    save('DGratings.mat','normF1','normF2','normSpikes')
    end
    end

elseif strcmp(params.stimName,'bars')
    disp('bars')

%This finds the position of bar intensity in the index holder to sort ON/OFF responses correctly for 
%radial plots
stringarrayforIntensity = string(splitCell(1,:));
intCode = stringarrayforIntensity=="intensity";
intCode = find(intCode > 0);

timings = timings*10;

figure(10);clf;
onResp =[];
offResp=[];
for a = 1:size(indexHolder,2)
        orientations = splitCell{2,size(splitCell,2)};
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
        clear onResp
        clear offResp
    
    for u = 1:length(unique(orientations))
        
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
            
            meanRate = mean(sum(spikeMatrix(sortedIndex(finalInd),1:timings(1)),2));
            
% meanRate=0;
         
            
            if double(indexHolder{1,a}(intCode)) > 0
                onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1):(timings(1)+timings(2))),2))-meanRate;
               
                offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2)):sum(timings)),2))-meanRate;
            else
                onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2)):sum(timings)),2))-meanRate;
          
%                 offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1):(timings(1)+timings(2))),2));
                offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1):(timings(1)+timings(2))),2))-meanRate;
            end
            
            %plot traces
            
            psthBGSubtracted = [];
            psthBGSubtracted = psthMatrix(sortedIndex(finalInd),:) - mean(psthMatrix(sortedIndex(finalInd),1:5000),2);
            psthData(u,:) = mean(psthBGSubtracted,1)-mean(psthBGSubtracted(1:5000),2);
           
%             psthData(u,psthData(u,:)<params.threshold)=0;
            
%             if double(indexHolder{1,a}(intCode)) > 0
%                 onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1):(timings(1)+timings(2))),2));
%                
%                 offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2)):sum(timings)),2));
%             else
%                 onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2)):sum(timings)),2));
%           
%                 offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1):(timings(1)+timings(2))),2));
%             end
            
%             psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(15); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
            
%             rawData(u,:) = mean(epochStorage(sortedIndex(finalInd),:),1);
%             figure(11); hold on
%             subplot(1,size(indexHolder,2),a)
%             plot(rawData(u,:)+(325*(u-1)))
%             title(indexHolder{1,a})
%             
            spikeData(u,:) = mean(spikeMatrix(sortedIndex(finalInd),:),1);
            figure(12); hold on
            subplot(1,size(indexHolder,2),a)
            plot(spikeData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
    end
if a == params.saveIter

    psthDataOut = psthData';
    psthDataOut = psthDataOut/max(max(psthDataOut)); %normalized!
save('oBarsPSTH.mat','psthDataOut')

end
% size(onResp)
% uniqueOrientations = unique(splitCell{2,size(splitCell,2)});
% uniqueOrientations
% onResp

if params.dataType == 1
    onResp(onResp<0)=0;
    offResp(offResp<0)=0;
end

if params.dataType < 1
    onResp = onResp/1000;
    offResp = offResp/1000;
end
uniqueOrientations = sort(str2double(uniqueOrientations));

allRadii = [uniqueOrientations (uniqueOrientations+180)];
allRadii = allRadii*2*pi/360;

figure(a)
subplot(1,2,1);
% polar(uniqueOrientations * 2 * pi /360, onResp);
polar(allRadii, [onResp onResp]);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
hold on
subplot(1,2,2);
polar(allRadii, [offResp offResp]);
title('OFF')
% xlabel(num2str(indexHolder{1,a}'));
xlabel(indexHolder{1,a}');

figure(a+size(indexHolder,2))
subplot(1,2,1);
polar(uniqueOrientations * 2 * pi /360, onResp);
% polar(allRadii, [onResp onResp]);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
hold on
subplot(1,2,2);
polar(uniqueOrientations * 2 * pi /360, offResp);
title('OFF')
% xlabel(num2str(indexHolder{1,a}'));
xlabel(indexHolder{1,a}');

if double(indexHolder{1,params.saveIter}(1)) <= 0 
    polarity = 0; % 1 on 0 off
else
    polarity = 1;
end
outOrientations = uniqueOrientations';
if a==params.saveIter
    
%     onRespOut = [onResp onResp]';
%     offRespOut = [offResp offResp]';
onRespOut=onResp';
offRespOut=offResp';
    save('oBarsResp.mat','onRespOut','offRespOut','polarity','outOrientations')
end

end
elseif strcmp(params.stimName,'moving bar')
    figure(10)
resp = [];
for a = 1:size(indexHolder,2)
        orientations = splitCell{2,size(splitCell,2)};
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
        clear resp
    
    for u = 1:length(unique(orientations))
        
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
         
            respBGSubtracted = [];
            respBGSubtracted = psthMatrix(sortedIndex(finalInd),params.recordTime) - mean(psthMatrix(sortedIndex(finalInd),1:5000),2);
%             respBGSubtracted(respBGSubtracted<0)=0;
            respBGSubtracted(respBGSubtracted<params.threshold)=0;
%             figure(11); clf;
%             plot(respBGSubtracted')
%             pause
%             hold on
%             plot(mean(respBGSubtracted,1))
%             pause
%             resp(u) = mean(mean(respBGSubtracted,2));
           
%             resp(u)=max(mean(respBGSubtracted,1));
            
%             
%             resp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),params.recordTime),2)) - mean(sum(spikeMatrix(sortedIndex(finalInd),1:5000),2));
%              resp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),params.recordTime),2));
            %plot traces
            
            psthBGSubtracted = [];
            psthBGSubtracted = psthMatrix(sortedIndex(finalInd),:) - mean(psthMatrix(sortedIndex(finalInd),1:5000),2);
            psthData(u,:) = mean(psthBGSubtracted,1);
            if params.dataType ==1
            psthData(u,psthData(u,:)<params.threshold)=0;
            psthScaler = 300;

            end
            resp(u)=mean(mean(psthData(u,params.recordTime),2));
            
%             psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
           
            
            if size(indexHolder,2) > 1 %subplot deletes first trace for some reason for current traces???  Or potentially when only 1 plot needed.  Don't understand, couldn't fix w/o if statement
            subplot(1,size(indexHolder,2),a)
            end
           
            plot(psthData(u,:)+(300*(u-1)))
            title(indexHolder{1,a})

%             rawData(u,:) = mean(epochStorage(sortedIndex(finalInd),:),1);
%             figure(11); hold on
%             subplot(1,size(indexHolder,2),a)
%             plot(rawData(u,:)+(325*(u-1)))
%             title(indexHolder{1,a})
%             
            spikeData(u,:) = mean(spikeMatrix(sortedIndex(finalInd),:),1);
            figure(14); hold on
            subplot(1,size(indexHolder,2),a)
            plot(spikeData(u,:)+(1*(u-1)))
            title(indexHolder{1,a})
%             
    end

if a==params.saveIter
    psthDataOut = psthData';
%     psthDataOut = psthDataOut/max(max(psthDataOut));
%     psthDataOut = psthDataOut/2;
    regionOut = zeros(size(psthMatrix,2),1);
    regionOut(params.recordTime(1):params.recordTime(size(params.recordTime,2)),1) = ones(length(params.recordTime),1);
        save('mBarPSTH.mat','psthDataOut','regionOut')
end
    
    

    
    
% if isequal(size(uniqueOrientations),size(resp))
%     disp('equal')
% holdOrientations = uniqueOrientations;
% elseif ~isequal(size(uniqueOrientations),size(resp))
% uniqueOrientations = holdOrientations;
% end


uniqueOrientations = sort(str2double(uniqueOrientations));

% allRadii = uniqueOrientations;
% allRadii = allRadii*2*pi/360;

size(uniqueOrientations);
size(resp);

resp

figure(a)
polar(uniqueOrientations * 2 * pi /360, resp);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
ylabel(indexHolder{1,a}');

figure(a+20)
bar(uniqueOrientations,resp)
outOrientations = uniqueOrientations';
if double(indexHolder{1,params.saveIter}(1)) <= 0 
    polarity = 0; % 1 on 0 off
else
    polarity = 1;
end

if a==params.saveIter
    respOut = resp';
    save('mBarsResp.mat','respOut','polarity','outOrientations')
end

end



end
end