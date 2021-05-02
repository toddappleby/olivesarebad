%popmikeOMD..it'sbad
cellName = 'A1'
typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == 'A');
    clear meanOut
    for z = 1:size(fileIndex,3)
    
    load(strtrim(fileIndex(:,:,z)))
    
       
        
        
        lister = unique(stimulusParams.spaceConstant);
        for y = 1:size(lister,2)
            
       tester(y) = mean(sum(analysisParams.response(stimulusParams.spaceConstant==lister(y),:),2));

        end
        tester=tester/max(tester);
        meanOut(z,:)=tester;
        
        
        
    end
    
    
        if size(meanOut,1)>1
            normAvg = mean(meanOut,1);
            normError = sem(meanOut,1);
        else
            normAvg = meanOut;
            normError = [];
        end
        
        outArray = normAvg';
        outError = normError';
        save(strcat('x',cellName),'outArray','outError')
            
    
        
        
        
    
    
    