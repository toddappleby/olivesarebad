function [importDataDir, exportForIgor] = getDirectory()

computerUser = getenv('USERNAME');

switch computerUser
    case 'reals'
        importDataDir = 'C:\Users\reals\Documents\PhD 2021\ClarinetExports\';
        exportForIgor = 'C:\Users\reals\Documents\PhD 2021\ClarinetExports\ImportToIgor';
    case 'todda'
        importDataDir = 'C:\Users\todda\Documents\Primate Data\ClarinetExports\';
        exportForIgor = 'C:\Users\todda\Documents\Primate Data\ImportToIgor\';
end