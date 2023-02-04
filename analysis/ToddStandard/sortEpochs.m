function sortedEpochs = sortEpochs(dataStruct)


for j = 1:size(dataStruct,2)
    
    str=string(dataStruct(j).meta.displayName);
    
    analysisArray(j,:) = str;
    %each epoch's analysis label 

    if j == 1
    holdUnique(j,:) = str;
    elseif isempty(find(holdUnique(:,:)==str,1))
    holdUnique(size(holdUnique,1)+1,:) = str;
    end
% holdUnique lists the analyses present  
end


% this gets the indices for each analysis
for i = 1:size(holdUnique,1)
    anArray = find(analysisArray == holdUnique(i));
    indexAnalysis(1:size(anArray,1),i) = anArray;
end

% this sets the analysis fields in the returned struct
for s = 1:size(indexAnalysis,2)
    spaceIt = holdUnique(s);
    spaceIt = char(spaceIt);
    spaceIt = spaceIt(~isspace(spaceIt));
    index = nonzeros(indexAnalysis(:,s));
    for h = 1:length(index)
    sortedEpochs(h).(spaceIt) = dataStruct(index(h)).epoch;
    end
end



    





