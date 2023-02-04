function [psthData] = orientedStim(epochStorage,spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

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
        orientations = str2double(splitCell{2,size(splitCell,2)});
        orientations = orientations((indexHolder{2,a}));
          [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
          sortedIndex = sortedIndex(I);
      

        for u = 1:length(unique(orientations))
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));

collectIndices = sortedIndex(finalInd);
            collectF1 = [];
            collectF2 = [];
            collectSpikes = [];
            holder = zeros(1,length(spikeMatrix));
            for uu = 1:length(sortedIndex(finalInd))
                holder(spikeMatrix(collectIndices(uu),:)>0)=sampleRate; %mike does this in his online analysis.  magnitude post binning is v small when this not done
                binnedData = binData(holder(stimStart:stimEnd),100,sampleRate);
                [F, phase] = frequencyModulation(binnedData, ...
                100, double(indexHolder{1,a}(1,tFreq)), 'avg', 1:2, []);
                collectF1 = [collectF1; F(1)]; 
                
                  
                 collectF2 = [collectF2; F(2)];
                
                
                collectSpikes = [collectSpikes; sum(spikeMatrix(collectIndices(uu),stimStart:stimEnd),2);];
            end
%
             

                avgF1(u) = mean(collectF1);
                avgF2(u) = mean(collectF2);
                spikeCount(u) = mean(collectSpikes,1);


%             data = mean(spikeMatrix(sortedIndex(finalInd),:),1);
%                     binnedData = BinSpikeRate(data(stimStart:stimEnd), 100, sampleRate);
%                 [F, phase] = frequencyModulation(binnedData, ...
%                 100, str2num(indexHolder{1,a}(1,tFreq)), 'avg', 1:2, params.killCycle1);
%                 avgF1(u)= F(1);
%                 avgF2(u)= F(2);     
%                 spikeCount(u) = mean(data(stimStart:stimEnd),2);
%                 
             psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
           
            figure(15); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(100*(u-1)))
            title(indexHolder{1,a})
        end

        
%         uniqueOrientations = sort(str2double(uniqueOrientations));
          
        uniqueOrientations = sort(uniqueOrientations);
       
        avgF1(1:length(uniqueOrientations));
    figure(a)
    subplot(1,3,1);
    polar(uniqueOrientations * 2 * pi/360, avgF1(1:length(uniqueOrientations)));
%     titleString = strcat({'F1 Resp'}, num2str(indexHolder{1,1}));
    title('F1 Resp')
    xlabel(analysisSplitters(1:3))
    hold on
    subplot(1,3,2);
    polar(uniqueOrientations * 2 * pi/360, avgF2(1:length(uniqueOrientations)));
%     titleString = strcat({'F2 Resp'}, num2str(indexHolder{1,1}));
    title('F2 Resp')
    xlabel(indexHolder{1,a}')
    
    subplot(1,3,3);
    polar(uniqueOrientations * 2 * pi/360, spikeCount(1:length(uniqueOrientations)));
%     titleString = strcat({'Spike Count'}, num2str(indexHolder{1,1}));
    title('Spike Count')
    xlabel(analysisSplitters(1:3))
    
    figure(10)

   
    
    end

elseif strcmp(params.stimName,'bars')
    disp('bars')

%This finds the position of bar intensity in the index holder to sort ON/OFF responses correctly for 
%radial plots
stringarrayforIntensity = string(splitCell(1,:));
intCode = stringarrayforIntensity=="intensity";
intCode = find(intCode > 0);

figure(10);clf;
onResp =[];
offResp=[];
for a = 1:size(indexHolder,2)
        orientations = str2double(splitCell{2,size(splitCell,2)});
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
        clear onResp
        clear offResp
    
    for u = 1:length(unique(orientations))
        
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
            
            meanRate = mean(sum(spikeMatrix(sortedIndex(finalInd),1:timings(1)*10),2));
         
            %take out mean rate by epoch
            
            for g = 1:length(finalInd)
                
                bgRate = sum(spikeMatrix(sortedIndex(finalInd(g)),1:timings(1)*10))/(timings(1)/10e2);
                onsetRate = sum(spikeMatrix(sortedIndex(finalInd(g)),timings(1)*10:(timings(1)+timings(2))*10))/(timings(2)/10e2);
                offsetRate = sum(spikeMatrix(sortedIndex(finalInd(g)),(timings(1)+timings(2))*10:sum(timings)*10))/(timings(3)/10e2);
%             offsetRate = sum(spikeMatrix(sortedIndex(finalInd(g)),(timings(1)+timings(2))*10:11500))/(timings(3)/10e2);


           
                subbedRateOn(g)= onsetRate - bgRate;
                subbedRateOff(g)= offsetRate - bgRate;
            end
                
            if params.dataType == 1
                if params.bgRate == 1
                    
                    
                    if double(indexHolder{1,a}(intCode)) > 0
                        onResp(u) = mean(subbedRateOn);
                        offResp(u) = mean(subbedRateOff);
                        
                        %                 onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
                        %                 offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));
                    else
                        onResp(u) = mean(subbedRateOff);
                        offResp(u) = mean(subbedRateOn);
                        
                        %                 onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));
                        %                 offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
                    end
                    
                else
                    
                    if double(indexHolder{1,a}(intCode)) > 0
                        %                 onResp(u) = mean(subbedRateOn);
                        %                 offResp(u) = mean(subbedRateOff);
                        
                        onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
                        offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));
                        
                    else
                        %                 onResp(u) = mean(subbedRateOff);
                        %                 offResp(u) = mean(subbedRateOn);
                        
                        onResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));
                        offResp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
                        
                    end
                end
            else
                
                if double(indexHolder{1,a}(intCode)) > 0
                    
                    onResp(u) = mean(mean(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
                    offResp(u) = mean(mean(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));
                    
                else
                    
                    onResp(u) = mean(mean(spikeMatrix(sortedIndex(finalInd),(timings(1)+timings(2))*10:sum(timings)*10),2));
                    offResp(u) = mean(mean(spikeMatrix(sortedIndex(finalInd),timings(1)*10:(timings(1)+timings(2))*10),2));
                    
                end
            end
           
               
               
            %plot traces
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(15); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(300*(u-1)))
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
            plot(spikeData(u,:)+(1*(u-1)))
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


% uniqueOrientations = sort(str2double(uniqueOrientations));
uniqueOrientations = sort(uniqueOrientations);

allRadii = [uniqueOrientations (uniqueOrientations+180)];
allRadii = allRadii*2*pi/360;
% allRadii = uniqueOrientations*2*pi/360;

figure(a)
subplot(1,2,1);
polar(uniqueOrientations * 2 * pi /360, onResp);
% % polar(allRadii, [onResp onResp]);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
hold on
subplot(1,2,2);
% polar(allRadii, [offResp offResp]);
polar(uniqueOrientations * 2 * pi /360, offResp);
title('OFF')
% xlabel(num2str(indexHolder{1,a}'));
xlabel(indexHolder{1,a}');

figure(a+20)
subplot(1,2,1);
polar([allRadii 0], [onResp onResp onResp(1)]);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
hold on
subplot(1,2,2);
polar([allRadii 0], [offResp offResp offResp(1)]);
title('OFF')
% xlabel(num2str(indexHolder{1,a}'));
xlabel(indexHolder{1,a}');

% figure(a+40)
% subplot(1,2,1);
% plot(uniqueOrientations, onResp)
% title('ON')
% xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
% hold on
% subplot(1,2,2);
% plot(uniqueOrientations, offResp)
% title('OFF')
% % xlabel(num2str(indexHolder{1,a}'));
% xlabel(indexHolder{1,a}');

figure(a+45)
plot(uniqueOrientations, onResp/max(abs(onResp)))
title('ON and OFF responses')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
hold on
plot(uniqueOrientations, offResp/max(abs(offResp)))
% xlabel(num2str(indexHolder{1,a}'));
xlabel(indexHolder{1,a}');

maxOn = max(abs(onResp));
maxOnOrientation = uniqueOrientations(onResp==maxOn);
minOffOpposite = offResp(onResp==maxOn);

maxOff = max(abs(offResp));
maxOffOrientation = uniqueOrientations(offResp==maxOff);
minOnOpposite = onResp(offResp==maxOff);
fileName = params.fileName;



if a==params.saveIter && params.saveGraph == 1
    if exist('C:\Users\reals\Documents\PhD 2021\ClarinetExports\PreferredOrientations.mat') == 2
    save('C:\Users\reals\Documents\PhD 2021\ClarinetExports\PreferredOrientations.mat','maxOn','maxOnOrientation','minOffOpposite','maxOff','minOnOpposite','fileName','-append')
    else
    save('C:\Users\reals\Documents\PhD 2021\ClarinetExports\PreferredOrientations.mat','maxOn','maxOnOrientation','minOffOpposite','maxOff','minOnOpposite','fileName')
    end
end

if double(indexHolder{1,params.saveIter}(1)) <= 0 
    polarity = 0; % 1 on 0 off
else
    polarity = 1;
end

if a==params.saveIter
    
    onRespOut = [onResp onResp]';
    offRespOut = [offResp offResp]';
    save('oBarsResp.mat','onRespOut','offRespOut','polarity')
end

end
elseif strcmp(params.stimName,'moving bar')
    figure(10);clf;
resp = [];
for a = 1:size(indexHolder,2)
        orientations = str2double(splitCell{2,size(splitCell,2)});
        orientations = orientations((indexHolder{2,a}));
         [orientations I] = sort(orientations);
         sortedIndex = indexHolder{2,a};
         sortedIndex = sortedIndex(I);
        clear resp
    
    for u = 1:length(unique(orientations))
        
            uniqueOrientations = unique(orientations); %because I don't know how to index the unique function
            finalInd = find(orientations==uniqueOrientations(u));
            resp(u) = mean(sum(spikeMatrix(sortedIndex(finalInd),params.Region),2));
            
            
            %plot traces
            psthData(u,:) = mean(psthMatrix(sortedIndex(finalInd),:),1);
            figure(10); hold on
            subplot(1,size(indexHolder,2),a)
            plot(psthData(u,:)+(600*(u-1)))
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
            plot(spikeData(u,:)+(1*(u-1)))
            title(indexHolder{1,a})
%             
    end

if a==params.saveIter
    psthDataOut = psthData';
    psthDataOut = psthDataOut/max(max(psthDataOut));
    regionOut = ones(length(params.Region),1);
        save('mBarPSTH.mat','psthDataOut','regionOut')
end
    
    

    
    
% if isequal(size(uniqueOrientations),size(resp))
%     disp('equal')
% holdOrientations = uniqueOrientations;
% elseif ~isequal(size(uniqueOrientations),size(resp))
% uniqueOrientations = holdOrientations;
% end


% uniqueOrientations = sort(str2double(uniqueOrientations));
uniqueOrientations = sort(uniqueOrientations);

% allRadii = uniqueOrientations;
% allRadii = allRadii*2*pi/360;

size(uniqueOrientations);
size(resp);



figure(a)
polar(uniqueOrientations * 2 * pi /360, resp);
title('ON')
xlabel(analysisSplitters(1:(length(analysisSplitters)-1)));
ylabel(indexHolder{1,a}');

% figure(a+20)
% bar(uniqueOrientations,resp)

if double(indexHolder{1,params.saveIter}(1)) <= 0 
    polarity = 0; % 1 on 0 off
else
    polarity = 1;
end

if a==params.saveIter
    respOut = resp';
    save('mBarsResp.mat','respOut','polarity')
end

end



end
end