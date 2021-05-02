cd('E:\Data Analysis_2020\2020_0420\')
load('20200420Bc2.mat')
%%

protocolToSort = 'Object Motion Texture';

counter = 0;
for z = 1:length(epochs)
    displayName = epochs(z).meta.displayName;
    if strcmp(displayName,protocolToSort)
        counter = counter +1;
        protocolIndices(counter) = z;
        epochTime(counter) = epochs(z).meta.epochTime;
        epochStartTime(counter) = epochs(z).meta.epochStartTime;
        stimulusClass(counter,:) = string(epochs(z).meta.stimulusClass);
    end
end

uniqueOrientations = unique(stimulusClass,'rows');

[yS,xI]= sortrows(stimulusClass);

figure(1)
scatter(sort(epochTime),categorical(yS))
% scatter(sort(epochTime),categorical(string(num2str(yS))))
figure(2)
scatter(sort(epochTime),xI)
% scatter(xI,categorical(yS))
% scatter(xI,categorical(string(num2str(yS))))

%% 

outEpochs = struct();
for c = 1:length(outData(:,1))
outEpochs(c).meta = epochs(protocolIndices(outData(c,2))).meta;
outEpochs(c).epoch = epochs(protocolIndices(outData(c,2))).epoch;
% outEpochs(c).uniqueCats = unique(uniqueOrientations(outData(c,2),:),'rows');
end

clear epochs
epochs = struct();
epochs = outEpochs;
save('Bc2_gratingDSOS_orientationlateBackwards.mat','epochs');

%%