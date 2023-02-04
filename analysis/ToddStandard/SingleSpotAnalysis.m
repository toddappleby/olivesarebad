function results = SingleSpotAnalysis(node,params)

%storage setup
temp = node{1}.parent.parent.parent.parent.splitValue;
index = find(temp == ' ');
titleString = temp(1:index);
titleString = strcat(strcat(titleString, 'TA-'), node{1}.parent.parent.parent.parent.splitValue);
title(titleString);
storageDir = strcat(params.rootDir, titleString);

if (~exist(storageDir, 'dir'))
     mkdir(storageDir);
end
cd(storageDir);

%timings setup
preTimeMS = node{1}.epochList.elements(1).protocolSettings.get('preTime');
preTime = preTimeMS / 1000; % now in seconds
stmTimeMS = node{1}.epochList.elements(1).protocolSettings.get('stimTime');
stmTime = stmTimeMS / 1000; %now in seconds
tailTimeMS = node{1}.epochList.elements(1).protocolSettings.get('tailTime');
tailTime = tailTimeMS / 1000; %now in seconds
samplingInterval = node{1}.epochList.elements(1).protocolSettings.get('sampleRate');
samplingInterval = 1/samplingInterval; %now in seconds/cycle
prePts = preTime / samplingInterval; %samples per length of time
stmPts = stmTime / samplingInterval;
tailPts = tailTime / samplingInterval;

allTime = preTime+stmTime+tailTime; %seconds
allPts = prePts + stmPts + tailPts;

sampleRate = node{1}.epochList.elements(1).protocolSettings.get('sampleRate');

for curNode = 1:length(node)

    
    nodeData = getSelectedData(node{curNode}.epochList,params.Amp);
    spikeMatrix = zeros(size(nodeData,1),size(nodeData,2));
   
    
    for curEpoch = 1:size(nodeData,1)
%           S = spikeDetectorOnline(wavefilter(nodeData(epoch,:),6));
%           spikeTimes(curEpoch) = S.sp;
          spikes = convertSpikes(nodeData(curEpoch,:));
          if isempty(spikes)
              disp('deleted epoch')
          else
          spikeMatrix(curEpoch,:) = spikes;
          psthMatrix(curEpoch,:) = psth(spikeMatrix(curEpoch,:),6+2/3,sampleRate,1);
          end
    end
    
%in case I deal with multiple nodes in the future
% nodeHolder = zeros(length(node),size(nodeData,2));
% size(spikeMatrix(curEpoch,:))
% nodeHolder(curEpoch) = spikeMatrix(curEpoch);
 
          results = spikeMatrix;
         %axis in seconds, but number of plotted points = samples. plot using prePts / stmPts as indices)
          xvals = linspace(0,allTime,allPts);
          
          figure(3); clf; hold on
          plot(xvals,mean(psthMatrix)) 
          lightON = zeros(1,size(xvals,2));
          lightON(1,prePts) = max(mean(psthMatrix)) + 10; 
          plot(xvals,lightON,'k-.')
              
          lightOFF = zeros(1,size(xvals,2));
          lightOFF(1,prePts+stmPts) = max(mean(psthMatrix)) + 10;       
          plot(xvals,lightOFF,'k-.')
%           
           axis([0 allTime 0 max(mean(psthMatrix))+12])
          title('Mean Spot Response from Darkness')
          ylabel('Spike Rate (Hz)')
          xlabel('time (sec)');
          
         figure(4);  clf; 
         plot(nodeData(1,:))
         title('example single spot epoch')
         xlabel('time points')
         ylabel('pA')
          
         if (params.SaveGraphs)
          
          if (params.CellAttached)
              figure(3)
          set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125],'Color','w')
         
          export_fig 'SSFig.pdf'
          figure(4)
          set(gcf,'Color','w')
          export_fig 'SSFig.pdf' -append
          else
              if node{1}.epochList.elements(1).protocolSettings.get('stimulus:Amp1:offset') < 1
              figure(4)
              set(gcf,'Color','w')
              export_fig 'SSFigWC.pdf'
              else
                  figure(4)
              set(gcf,'Color','w')
              export_fig 'SSFigWC.pdf' -append
          end
          
          end
         end
end

