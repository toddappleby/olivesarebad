%% 
% exportFolder = '/Users/toddappleby/Documents/Data/Clarinet Exports/';
exportFolder = getDirectory();
protocolToAnalyze = 'SaccadeAndPursuit';
%typesFolder ='C:\Users\reals\Documents\PhD 2023\ClarinetExports\MikeSaccade'
typesFolder = 'C:\Users\todda\Documents\Primate Data\ClarinetExports\MikeSaccade'

cd(typesFolder)
dataDir = dir;
%figure out which folders start with 2! (data folders), then create array of
%folder names
folderSet = char(string({dataDir.name}));
% folderSet(:,:,1:3)=[]; %removes macOS hidden stuff
folderSet(:,:,1:2)=[]; 

%%     
    ParaOrSmooth = 'S';

    typeDir = dir; 
    fileSet = char(string({typeDir.name})); 
    fileIndex = fileSet(:,:,fileSet(:,1,:) == ParaOrSmooth);

    popCRSacc=zeros(1,13);
         popCRStatic=zeros(1,13);
                popCRErrSacc=zeros(1,13);
                popCRErrStatic=zeros(1,13);
    
    for g =1:size(fileIndex,3)
%         
%      
        load(fileIndex(:,:,g))

        
        if strcmp(ParaOrSmooth,'P') && strcmp(fileIndex(:,:,g),'ParasolOn_crf_extracellular_20180125Ac1_grating.mat')
        zeroContrastPSTHsacc(g,1:size(avgData,3)) = avgData(2,1,:);
        zeroContrastPSTHstatic(g,1:size(avgData,3)) = avgData(3,1,:);    
        else
        zeroContrastPSTHsacc(g,1:size(avgData,3)) = avgData(1,7,:);
        zeroContrastPSTHstatic(g,1:size(avgData,3)) = avgData(2,7,:);
        end
        
        if isfield(analysisParams,'resp')
        contrastResponse = analysisParams.resp;
%         contrastResponse(1,:) = contrastResponse(1,:)-contrastResponse(1,7);
%         contrastResponse(2,:) = contrastResponse(2,:)-contrastResponse(2,7);
%         contrastResponse(contrastResponse<0)=0;
        contrastResponse(1,:) = contrastResponse(1,:)/max(contrastResponse(1,:));
        contrastResponse(2,:) = contrastResponse(2,:)/max(contrastResponse(2,:));
        contrastResponseErr = analysisParams.err;
        contrastRespX = stimulusParams.uniqueCt;
        else
        contrastResponse = resp(1:2,:);
%         contrastResponse(1,:) = contrastResponse(1,:)-contrastResponse(1,7);
%         contrastResponse(2,:) = contrastResponse(2,:)-contrastResponse(2,7);
%         contrastResponse(contrastResponse<0)=0;
        contrastResponse(1,:) = contrastResponse(1,:)/max(contrastResponse(1,:));
        contrastResponse(2,:) = contrastResponse(2,:)/max(contrastResponse(2,:));
        contrastResponseErr = zeros(2,7);
        contrastRespX = stimulusParams.uniqueCt;
%             contrastResponse(1:2,7:13) = resp(2:3,:);
%             contrastResponseErr(1:2,7:13)=ones(2,size(resp,2));
%            contrastRespX(7:13) = stimulusParams.uniqueCt;
%            contrastRespX(1:6)= -[1,.75,.5,.25,.125,.0625];
        end
%       errorbar([contrastRespX;contrastRespX],contrastResponse,contrastResponseErr)
        
%         for z = 1:size(avgData,2)
%             holdover = avgData(1,z,:);
%             averageDataSacc(z,:) = holdover(:)';
%             holdover2 = avgData(2,z,:);
%             averageDataStatic(z,:) = holdover2(:)';
%         end

%         figure(3)
% errorbar(contrastRespX,contrastResponse(1,:),contrastResponseErr(1,:))
% hold on
% errorbar(contrastRespX,contrastResponse(2,:),contrastResponseErr(2,:))

        if length(stimulusParams.uniqueCt) == 7 && strcmp(ParaOrSmooth,'S')
                popCRSacc(g,7:13) = contrastResponse(1,:); 
                popCRStatic(g,7:13) = contrastResponse(2,:);
                popCRErrSacc(g,7:13) = contrastResponseErr(1,:);
                popCRErrStatic(g,7:13) = contrastResponseErr(2,:);
        elseif length(stimulusParams.uniqueCt) == 7 && strcmp(ParaOrSmooth,'P')
                popCRSacc(g,4:10) = contrastResponse(1,:); 
                popCRStatic(g,4:10) = contrastResponse(2,:);
                popCRErrSacc(g,4:10) = contrastResponseErr(1,:);
                popCRErrStatic(g,4:10) = contrastResponseErr(2,:);
        else
                popCRSacc(g,:) = contrastResponse(1,:); 
                popCRStatic(g,:) = contrastResponse(2,:);
                popCRErrSacc(g,:) = contrastResponseErr(1,:);
                popCRErrStatic(g,:) = contrastResponseErr(2,:);
         end
 end
    
    
subplot(2,1,1);plot(zeroContrastPSTHsacc')
subplot(2,1,2);plot(zeroContrastPSTHstatic')
figure;plot(reshape(avgData(1,:,:),13,size(avgData(1,7,:),3))')
%% manual load and plot
% numLines = 1;
% lineToUse = 7;
% plot(reshape(avgData(1,lineToUse,:),numLines,size(avgData(1,lineToUse,:),3))');hold on; plot(reshape(avgData(2,lineToUse,:),numLines,size(avgData(1,lineToUse,:),3))');line([1250 1250],[0 130],'Color','k','LineStyle','--');line([analysisParams.samplePts(1) analysisParams.samplePts(1)],[0 130],'Color','r','LineStyle','--');line([analysisParams.samplePts(2) analysisParams.samplePts(2)],[0 130],'Color','r','LineStyle','--');legend('moving surround','static surround','motion end','sample range')
curDir = dir;
fileToLoad = curDir(3).name;
plotLoad = load(fileToLoad);
load(fileToLoad)

if isfield(plotLoad,'blockParams')
    moveEnd = plotLoad.blockParams.preTime + plotLoad.blockParams.waitTime;
    moveEnd = moveEnd;
else
    moveEnd = 1250;
end
    

all_contrasts = stimulusParams.uniqueCt;

pos_contrasts = all_contrasts(all_contrasts>0);

zeroIndex = find(all_contrasts==0);



numLines = 1;
lineToUse = 8:13;
count = 0;
for xx = zeroIndex:zeroIndex+length(pos_contrasts)
    count =count+1;
subplot(length(pos_contrasts)+1,1,count)
plot(reshape(avgData(1,xx,:),numLines,size(avgData(1,xx,:),3))');hold on; plot(reshape(avgData(2,xx,:),numLines,size(avgData(1,xx,:),3))');line([moveEnd moveEnd],[0 200],'Color','k','LineStyle','--','LineWidth',2);line([analysisParams.samplePts(1) analysisParams.samplePts(1)],[0 200],'Color','r','LineStyle','--','LineWidth',2);line([analysisParams.samplePts(2) analysisParams.samplePts(2)],[0 200],'Color','r','LineStyle','--','LineWidth',2);
if xx == 7
    legend('moving surround','static surround','motion end','sample range','FontSize',12)
end
title(all_contrasts(xx))
end
sgtitle('On Smooth Motion & Static Surround Positive Contrast Example Cell 1')
%% rename
curDir = dir;
replacementSampleNums = [7600 8350;5100 6000;5100 6000;1300 1460;1300 1460];
for p = 3:length(curDir)
    load(curDir(p).name)
    analysisParams.samplePts(1) = replacementSampleNums(p-2,1)
    analysisParams.samplePts(2) = replacementSampleNums(p-2,2)
    save(curDir(p).name,'analysisParams','-append')
end


%% :) 
% zeroContrastCondition(1,:) = popCRSacc(:,7);
% zeroContrastCondition(2,:) = popCRStatic(:,7);
popCRSacc(popCRSacc==0) = NaN;
popCRStatic(popCRStatic==0) = NaN;

saccMean = mean(popCRSacc,'omitnan');
staticMean = mean(popCRStatic,'omitnan');

saccMean(isnan(saccMean))=0;
staticMean(isnan(staticMean))=0;

saccErr = sem(popCRSacc,'omitnan');
staticErr = sem(popCRStatic,'omitnan');

saccErr(isnan(saccErr))=0;
staticErr(isnan(staticErr))=0;

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

%% hz shift vs gain shift in CRF curves

popCRSacc(isnan(popCRSacc))=0;
popCRStatic(isnan(popCRStatic))=0;
contrastRespX = contrastRespX';

for pp = 1:size(popCRSacc,1)
    paramsHz(pp,:) = fitMultiVarParams([contrastRespX;contrastRespX],[popCRSacc(pp,:);popCRStatic(pp,:)],1,[0,-1.5])
    modeledCRFHz = multiHZNL(paramsHz(pp,:),[contrastRespX;contrastRespX]);
%     figure;clf;plot(contrastRespX,modeledCRFHz(1,:),'r');hold on
%     plot(contrastRespX,modeledCRFHz(2,:),'k')
%     plot(contrastRespX,popCRSacc(pp,:),'r--')
%      plot(contrastRespX,popCRStatic(pp,:),'k--')
%      

    paramsGain(pp,:) = fitMultiVarParams([contrastRespX;contrastRespX],[popCRSacc(pp,:);popCRStatic(pp,:)],0,[0,-1.5]);
    modeledCRFGain = multiGainNL(paramsGain(pp,1:4),[contrastRespX;contrastRespX]);
%     figure(10);clf;plot(contrastRespX,modeledCRFGain(1,:),'r');hold on
%     plot(contrastRespX,modeledCRFGain(2,:),'k')
%      plot(contrastRespX,popCRSacc(pp,:),'r--')
%      plot(contrastRespX,popCRStatic(pp,:),'k--')
%     

errSaccHz(pp) = immse(popCRSacc(pp,:),modeledCRFHz(1,:));
errSaccGain(pp) = immse(popCRSacc(pp,:),modeledCRFGain(1,:));
errStaticHz(pp) = immse(popCRStatic(pp,:),modeledCRFHz(2,:));
errStaticGain(pp) = immse(popCRStatic(pp,:),modeledCRFGain(2,:));

hzParam(pp,:) = [paramsHz(pp,3),paramsHz(pp,4)];
gainParam(pp,:) = [paramsGain(pp,3),paramsGain(pp,4)];
    
end

hz1=(diff(hzParam(:,[1,2]),1,2)./sum(hzParam(:,[1,2]),2))*100;
gain1=(diff(gainParam(:,[1,2]),1,2)./sum(gainParam(:,[1,2]),2))*100;

plot(ones(1,size(hz1,1)),hz1,'.'); hold on
plot(ones(1,size(hz1,1))+.2,gain1*-1,'.')
axis([.5 1.5 -20 20])

plot(ones(1,size(hz1,1)),errSaccHz,'r.','MarkerSize',15); hold on
plot(1,mean(errSaccHz),'r.','MarkerSize',40)
plot(ones(1,size(hz1,1))+.2,errSaccGain,'g.','MarkerSize',15)
plot(1.2,mean(errSaccGain),'g.','MarkerSize',40)
errorbar(1,mean(errSaccHz),sem(errSaccHz),'r.','MarkerSize',40)
errorbar(1.2,mean(errSaccGain),sem(errSaccGain),'g.','MarkerSize',40)

axis([.5 1.5 0 .2])


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