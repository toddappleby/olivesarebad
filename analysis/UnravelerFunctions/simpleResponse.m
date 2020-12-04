function simpleResponse(spikeMatrix,psthMatrix,timings,splitCell,indexHolder,params)

sampleRate = 10000;
stimStart = (timings(1)*1e-3)*sampleRate+1;
stimEnd = (timings(1) + timings(2))*1e-3*sampleRate;
stimOff = (timings(1)+timings(2)+timings(3))*10;

figure(11)
subplot(2,1,1)

xvals = linspace(0,stimOff/10,length(psthMatrix));



for g = 1:size(indexHolder,2)

subplot(length(indexHolder),1,g)
plot(xvals,mean(psthMatrix(indexHolder{2,g},:),1))
xline(stimStart/10,'LineStyle','--','Color','r','LineWidth',2)
xline(stimEnd/10,'LineStyle','--','Color','r','LineWidth',2)
ylabel(indexHolder{1,g})
xlabel('time (ms)')

end

if strcmp(params.protocolID,'Single Spot')

spotThing = mean(psthMatrix(indexHolder{2,1},:),1);
spotThing = spotThing';
% saveName = strcat('spotResponse_',params.cellName);
save('spotResponse.mat','spotThing')
elseif strcmp(params.protocolID,'Led Pulse')
    
uvResponse = mean(psthMatrix(indexHolder{2,1},:),1);
uvResponse = uvResponse';

redResponse = mean(psthMatrix(indexHolder{2,2},:),1);
redResponse = redResponse';   
save('ledResponse.mat','uvResponse','redResponse')
else
    
uvResponse = mean(psthMatrix(indexHolder{2,1},:),1);
uvResponse = uvResponse';

redResponse = mean(psthMatrix(indexHolder{2,2},:),1);
redResponse = redResponse';    
save('chromaticResponse.mat','uvResponse','redResponse')
end
    
for t = 1:length(indexHolder)
    

figure(20)   
% spikedThing = sum(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1)/max(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1)),2);
spikedThing = sum(mean(spikeMatrix(indexHolder{2,t},params.moveTime:end),1));
scatter(t,spikedThing,50,'filled')
xticks([1 2 3 4])
labelHolder(1,t)=indexHolder{1,t};
hold on

end
xticklabels(labelHolder);



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