
for a = 1:length(Data)
    protocolList(a) = string(Data(a).protocols{2,2});
    labels(a) = string(Data(a).displayProperties{10,2});
    
end

availableProtocols = unique(protocolList);
cellList = unique(labels);
%Data(1).rigConfig{1,5}{2}{11}  -- ndf 
%Data(1500).rigConfig{5}{2}{3} -- offset
%find(test=='.')

%cell2struct(metaData(:,2:end),string(metaData(:,1)),1) make structs 
%%
for b = 1:length(cellList)
    
    %grab by cell
    currCellSelect = labels == cellList(b);
    
        epochs = struct();
        newData = Data(currCellSelect);
    for c = 1:length(newData)
        
        epochs(c).epoch = newData(c).responses.y{1,1}';
        
        
        
%         if(size(newData(c).responses.y,2)==1)
%      
%             continue
%             
%             if c == length(newData) 
%                 break
%             end
%             c
%         end
%         epochs(c).epoch = newData(c).responses.y{1,2}'; 
      
        

        metaData = cell(0,0);
        metaData = newData(c).protocols(:,1); 
        metaData = [metaData newData(c).protocols(:,2)];
        
        thingSize = (size(metaData,1)+size(newData(c).rigConfig{1},1));
        
        metaData(1:thingSize,1) = [metaData(:,1); newData(c).rigConfig{1}] ;
        
        metaData(1:thingSize,2) = [metaData(1:(thingSize-size(newData(c).rigConfig{2},1)),2); newData(c).rigConfig{2}];
        
        metaData(size(metaData,1)+1,1) = {'NDF'};
        metaData(size(metaData,1),2) = {newData(c).rigConfig{1,5}{2}{11}};
        
        metaData(size(metaData,1)+1,1) = {'centerOffsetTrue'};
        
        
        cOffset = newData(c).rigConfig{1,5}{2}{3}';
        cOffset = sprintf('%d%d%d', cOffset(1),cOffset(2));
        
        metaData(size(metaData,1),2) = {cOffset};
        
        metaData(2,1) = {'displayName'};
        
        getProtocolName = char(metaData(2,2));
        proName=getProtocolName(max(find(getProtocolName=='.')+1):end);
        metaData(2,2) = {proName};
     
        
        epochs(c).meta = cell2struct(metaData(:,2),string(metaData(:,1)),1);



    
      

    end
    save(strcat(cellList(b),'.mat'),'epochs')
end

%% original nori unraveler -- better but don't have time to use

for b = 1:length(cellList)
    
    %grab by cell
    currCellSelect = labels == cellList(b);
    

    for c = 1:length(availableProtocols)

        %grab by protocl
        logicalSelect = protocolList == availableProtocols(c);
        
        %combine (note: can't find epoch group labels.... hope can compensate by storing NDF and offset)
        finalProtocol = currCellSelect./logicalSelect;
        finalProtocol = finalProtocol == 1;
        
        if isempty(find(finalProtocol>0))
            continue
        end
        
        newData = Data(finalProtocol);
        sizeCheck = 0;
        for d = 1:length(newData)
            sizeCheck(d) = length(newData(1).responses.x{1,1});
        end
        xData = zeros(d,max(sizeCheck));
        yData = zeros(d,max(sizeCheck));
        xFrame = zeros(d,max(sizeCheck));
        yFrame = zeros(d,max(sizeCheck));
        metaData = cell(0,0);
        metaData = newData(1).protocols(:,1);


            for e = 1:length(newData) 
                
                xData(e,:) = newData(e).responses.x{1,1};
                yData(e,:) = newData(e).responses.y{1,1};
                if size(newData(e).responses.x,2) > 1
                xFrame(e,:) = newData(e).responses.x{1,2};
                yFrame(e,:) = newData(e).responses.y{1,2};
                end
                NDFstorage(e) = newData(e).rigConfig{1,5}{2}{11};
                centerOffset(e,:) = newData(e).rigConfig{5}{2}{3}';
                metaData = [metaData newData(e).protocols(:,2)];
                        
            end

      getProtocolName = char(availableProtocols(c));
      proName=string(getProtocolName(max(find(getProtocolName=='.')+1):end));
            
      if size(newData(e).responses.x,2) > 1  
      save(strcat(cellList(b),'_',proName,'.mat'),'xData','yData','xFrame','yFrame','metaData','centerOffset','NDFstorage')
      else
      save(strcat(cellList(b),'_',proName,'.mat'),'xData','yData','metaData','centerOffset','NDFstorage')
      end
    end
end

