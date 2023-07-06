%popmikeOMGratings (simple)

dataDir=dir;
folderSet = char(string({dataDir.name}));
folderSet(:,:,1:2)=[];                            

saveDir = 'C:\Users\todda\Documents\Primate Data\Saves\';

for c = 1:size(folderSet,3)

cellName = folderSet(:,:,c);  

cd(strcat('C:\Users\todda\Documents\Primate Data\Mike Data\ObjectMotionTexture\',folderSet(:,:,c))) 

typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileSet(:,:,1:2)=[];
%     fileIndex = fileSet(:,:,fileSet(:,1,:) == 'P' | 'S');
    fileIndex = fileSet;
    
outEyeObj = [];
outEye = [];
outObject = [];
outDrift = [];

outEyeObjMean = [];
outEyeMean = [];
outObjectMean = [];
outDriftMean = [];

    for z = 1:size(fileIndex,3)
    
    load(strtrim(fileIndex(:,:,z)))
    
        for y = 1:4
            meanCell(y) = mean(avgData(y,:),2);
            if exist('epochParams')
             classLabel(y) = string(epochParams(y).stimulusClass);
            end
        end

    outEyeObj(z) = meanCell(1);;
    outEye(z) = meanCell(2);;
    outObject(z) = meanCell(3);;
    outDrift(z) = meanCell(4);;


    outEyeObjMean(z) = meanCell(1);
    outEyeMean(z) = meanCell(2);
    outObjectMean(z) = meanCell(3);
    outDriftMean(z) = meanCell(4); 

    
 
    end

  differentialIndex{1,c}= (outEyeObj - outObject )./(outEyeObj + outObject);
  synchIndex{1,c} = (outEye - outObject)./(outEye + outObject);

%   differentialIndex(c,:) = (outEyeObj - outObject )./(outEyeObj + outObject);
%   synchIndex(c,:) = (outEye - outObject)./(outEye + outObject);
  cellList(c) = string(cellName);


  

end

save(strcat(saveDir,'OMG_',cellName),'differentialIndex','synchIndex','cellList','classLabel')

% 
colorWheel = ['k','g','b','r','cy']
% for g = 1:4 %always 4 conditions 
%     for p = 1:length(differentialIndex{1,1})
% 
% 
%         plot([ones(size(differentialIndex{1,g},2),1)]',differentialIndex{1,g},'.','MarkerSize',20,'Color',colorWheel(g))
%         hold on
%         plot(1,mean(differentialIndex{1,g}),'.','MarkerSize',30,'Color',colorWheel(g))
%     
%         plot([ones(size(differentialIndex{1,g},2),1)*2]',synchIndex{1,g},'.','MarkerSize',20,'Color',colorWheel(g))
%         hold on
%         plot(2,mean(synchIndex{1,g}),'.','MarkerSize',30,'Color',colorWheel(g))
%         
%     
%      end
% end

line([1 2],[0 0],'LineStyle','--','Color','k')

figure(98)
for g = 1:size(differentialIndex,2)
%     for p = 1:length(differentialIndex{1,1})


        plot(1,mean(differentialIndex{1,g}),'.','MarkerSize',30,'Color',colorWheel(g))

    hold on
    
        plot(1,sem(differentialIndex{1,g}),'.','MarkerSize',30,'Color',colorWheel(g))


        plot(2,mean(synchIndex{1,g}),'.','MarkerSize',30,'Color',colorWheel(g))

        plot(2,sem(synchIndex{1,g}),'.','MarkerSize',30,'Color',colorWheel(g))

        
    
%      end
end

line([1 2],[0 0],'LineStyle','--','Color','k')




% plot([1 1 1],differentialIndex,'.','MarkerSize',20)
% hold on
% plot([1 1 1],mean(differentialIndex,2),'.','MarkerSize',40)
% plot([2 2 2],synchIndex,'.','MarkerSize',20)
% plot([2 2 2],mean(synchIndex,2),'.','MarkerSize',40)
% 
% axis([0 3 -15 15])

  %% old 
% outEyeObj = outEyeObj/max(outEyeObj);
% outEye = outEye/max(outEye);
% outObject = outObject/max(outObject);

% outArray = [mean(outEyeObj/max(outEyeObj)) mean(outEye/max(outEye)) mean(outObject/max(outObject))]';
% outError = [sem(outEyeObj/max(outEyeObj)) sem(outEye/max(outEye)) sem(outObject/max(outObject))]';




% plot([mean(outEyeObj) mean(outEye) mean(outObject)])
% figure
% plot([sem(outEyeObj) sem(outEye) sem(outObject)])

% plot([outEyeObj outEye outObject])
% figure
% plot([outEyeObj outEye outObject])

plot(outArray)

outArray = flip(outArray);
outArray = outArray/max(outArray);

outError = flip(outError);
% outError = outError/max(outError);



% save(strcat('omG_',cellName),'outArray','outError')
    
        
        
        
    
    
    