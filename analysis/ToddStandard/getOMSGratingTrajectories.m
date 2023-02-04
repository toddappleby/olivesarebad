%Create trajectory for OMS grating

jitterSpeed = 1000; %pixel/second -- what is this in micron?  probably close to 1 wtf? doesn't make sense!
frameRate = 60;

secperFrame = 1/60 * 1000;

stepSize = 1000/60; 

jitterPerFrame = jitterSpeed*secperFrame;

barWidth = 100; %microns...100 is what I usually use.  barely matters for this

%seed it
%% grab params/info from data unraveler 

for z = 1:length(surroundSeed)
    
    stimClass(z) = splitCell{2,1}(z);
    
end
    
allClasses = unique(stimClass);
for y = 1:length(unique(stimClass))
    
    seedSet = surroundSeed(strcmp(allClasses(y),stimClass));
    
    if strcmp(allClasses(y),"eye+object")
       
        secondSeed = seedSet-1781; %because seed is surround seed which is +1781 in stim code
        
    end
    
    %get trajectories 
   
    
    for x = 1:length(seedSet)
        
        noiseStream = RandStream('mt19937ar', 'Seed', seedSet(x));
        
        for aa = 1:300 %num frames given 5000 ms stim time
            fullTraj(aa) = noiseStream.randn*2*pi * stepSize / 100;
        end
        
        if strcmp(allClasses(y),"eye+object")
       
         noiseStream2 = RandStream('mt19937ar', 'Seed', secondSeed(x));
         
         for bb = 1:300 %num frames given 5000 ms stim time
            fullTrajSurround(bb) = noiseStream2.randn*2*pi * stepSize / 100;
         end
       
         
        end
    
    
    
    end
end
