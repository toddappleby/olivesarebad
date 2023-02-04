expDate = '2019_1107';
cellNum = 'Bc1';


dataDir = 'C:\Users\reals\Documents\PhD 2021\ClarinetExports\';

cd(strcat(dataDir,expDate))
cellDate = strcat(expDate(1:4),expDate(6:9));

load(strcat(cellDate,cellNum))
cellData = epochs;
load(strcat(cellDate,cellNum,'_FT'))
frameTs = epochs;

for z = 1:size(epochs,2)
   list(z) = string(epochs(z).meta.displayName);
end
uniqueProtocols = unique(list')
%%

for i = 1:length(uniqueProtocols)
    
    epochs( 


