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
    fileIndex = fileSet(:,:,fileSet(:,1,:) == 'S');
    
    for g = 1:size(fileIndex,3)
        load(fileIndex(:,:,g))
           
        contrastResponseErr = analysisParams.err;
        contrastRespX = stimulusParams.uniqueCt;
        contrastResponse = analysisParams.resp;
%       errorbar([contrastRespX;contrastRespX],contrastResponse,contrastResponseErr)
        
        for z = 1:size(avgData,2)
            holdover = avgData(1,z,:);
            averageDataSacc(z,:) = holdover(:)';
            holdover2 = avgData(2,z,:);
            averageDataStatic(z,:) = holdover2(:)';
        end

%         figure(3)
% errorbar(contrastRespX,contrastResponse(1,:),contrastResponseErr(1,:))
% hold on
% errorbar(contrastRespX,contrastResponse(2,:),contrastResponseErr(2,:))


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
errorbar(contrastRespX,saccMean,saccErr);hold on
errorbar(contrastRespX,staticMean,staticErr)


