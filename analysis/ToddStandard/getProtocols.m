function uniqueProtocols = getProtocols(epochs)

uniqueProtocols = [];

 for z = 1:length(epochs)
     uniqueProtocols = [uniqueProtocols string(epochs(z).meta.displayName)];
 end
 uniqueProtocols = unique(uniqueProtocols);
end
% 
% uniqueProtocols = [];
% 
% for z = 1:size(epochs,2)
%    list(z) = string(epochs(z).meta.displayName);
%    
% %    allEpochData(z,1:length(epochs(z).epoch)) = epochs(z).epoch;
% end
% %find 0s, make index, create new string
% list;
% while ~isempty(list)
% uniqueProtocols = [uniqueProtocols; list(1)];%OK<AGROW>
% uniqueCheck = strcmp(uniqueProtocols(length(uniqueProtocols)),list);
% newIndex = find(uniqueCheck==0);
% list = list(newIndex);
% currentIndex = find(uniqueCheck==1);
% rawData{length(uniqueProtocols),1} = uniqueProtocols(length(uniqueProtocols));
% end
% 
% uniqueProtocols = sort(uniqueProtocols); %#ok<NOPTS>