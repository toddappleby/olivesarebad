%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = '/Users/reals/Documents/PhD 2021/ClarinetExports/';
protocolToAnalyze = 'SaccadeAndPursuit';
typesFolder ='C:\Users\reals\Documents\PhD 2021\ClarinetExports\MikeSaccade';

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
% folderSet(:,:,1:3)=[]; %removes macOS hidden stuff
folderSet(:,:,1:2)=[]; 

%%     
    

    typeDir = dir; 
    fileSet = char(string({typeDir.name})); 
    fileIndex = fileSet(:,:,fileSet(:,1,:) == 'P');
    
    for g = 1:size(fileIndex,3)
        load(fileIndex(:,:,g))
           
        contrastResponseErr = analysisParams.err;
        contrastRespX = stimulusParams.uniqueCt;
%         contrastResp = 
        if isfield(analysisParams,resp)
        contrastResponse = analysisParams.resp;
        else
            contrastResponse = resp;
        end
            
%       errorbar([contrastRespX;contrastRespX],contrastResponse,contrastResponseErr)
        
%         for z = 1:size(avgData,2)
%             holdover = avgData(1,z,:);
%             averageDataSacc(z,:) = holdover(:)';
%             holdover2 = avgData(2,z,:);
%             averageDataStatic(z,:) = holdover2(:)';
%         end

%         figure(3)
errorbar(contrastRespX,contrastResponse(1,:),contrastResponseErr(1,:))
hold on
errorbar(contrastRespX,contrastResponse(2,:),contrastResponseErr(2,:))
pause 

        popCRSacc(g,:) = contrastResponse(1,:);
        popCRStatic(g,:) = contrastResponse(2,:);
        popCRErrSacc(g,:) = contrastResponseErr(1,:);
        popCRErrStatic(g,:) = contrastResponseErr(2,:);
        
    end
    
%% :) 

saccMean = mean(popCRSacc);
staticMean = mean(popCRStatic);

saccErr = sem(popCRSacc);
staticErr = sem(popCRStatic);

figure
errorbar(contrastRespX,saccMean,saccErr,'r');hold on
errorbar(contrastRespX,staticMean,staticErr,'b')
% 
% figure(10)
% plot(stimulusParams.uniqueCt,popCRSacc(4,:));hold on
% plot(stimulusParams.uniqueCt,popCRStatic(4,:));

figure(10)
plot(stimulusParams.uniqueCt,popCRSacc(2,:));hold on
plot(stimulusParams.uniqueCt,popCRStatic(2,:));
% figure(11)


%% sensitivity measures

% jenson-shannon 
collectionJS=[];
collectionDSacc=[];
collectionDStatic = [];
contrastValues = stimulusParams.uniqueCt;
for c = 7:size(popCRSacc,2)
posContrastSacc = popCRSacc(:,c);
posContrastStatic = popCRStatic(:,c);
popCountCRSacc = posContrastSacc.*.1;
popCountCRStatic = posContrastStatic.*.1;
[Z, X]=analysis.stats.probabilitiesFromCounts(popCountCRSacc(:),popCountCRStatic(:));

figure(55);clf
plot(Z)
hold on
plot(X)
legend('sacc','static')
pause

JSresult = analysis.stats.jsDivergence(Z,X);

collectionJS(c-6) = JSresult;
end
figure(2)
plot(contrastValues(7:end),collectionJS)

%discriminability 

zeroContrastSacc = popCRSacc(:,7); % for full range of contrasts -- pick new 0 index for parasols with smaller contrast range
zeroContrastStatic = popCRStatic(:,7);

for c = 8:size(popCRSacc,2)
posContrastSacc = popCRSacc(:,c);
posContrastStatic = popCRStatic(:,c);
% popCountCRSacc = posContrastSacc.*.1;
% popCountCRStatic = posContrastStatic.*.1;
% dPrimeSacc=[];
% dPrimeStatic=[];
% for d = 1:length(posContrastSacc)
dPrimeSacc = analysis.stats.dprime(posContrastSacc,[zeroContrastSacc zeroContrastStatic]);
dPrimeStatic = analysis.stats.dprime(posContrastStatic,[zeroContrastSacc zeroContrastStatic]);
% end

collectionDSacc(c-7) = dPrimeSacc;
collectionDSstatic(c-7) = dPrimeStatic;
end
figure(4)
plot(contrastValues(8:end), collectionDSacc)
hold on
plot(contrastValues(8:end), collectionDSstatic)


% dPrimeSacc = analysis.stats.dprime(posContrastSacc(d),[zeroContrastSacc(d) zeroContrastStatic]);
% dPrimeStatic = analysis.stats.dprime(posContrastStatic(d),[zeroContrastSacc zeroContrastStatic(d)]);