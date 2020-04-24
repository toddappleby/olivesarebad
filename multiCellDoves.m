%% this is for loading images

cellName = '20200420Bc1_ONOFF';

load(['DovesData',cellName,'_Doves.mat'])

% for k=1:length(stimIndex)

fileId = fopen(['E:\Data Analysis_2020\code\Manookin Repository\manookin-package\resources\doves\images\',stimIndex{2,k}],'rb','ieee-be');
img = fread(fileId, [1536 1024], 'uint16');
fclose(fileId);
img = double(img');
img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level
backgroundIntensity = mean(img(:));%set the mean to the mean over the image
img = img.*255; %rescale s.t. brightest point is maximum monitor level
imageMatrix = uint8(img);

subplot(2,2,1:2)
% set(gcf,'position',[0,0,1280,920])
% axes('Position', [.5 1 .5 .5]);
colormap(gray(256));
imagesc(imageMatrix);
set(gca,'visible','off')
subplot(2,2,3)  
% axes('Position', [0.1 0.2 0.3 0.5]);
plot(stimIndex{3,k})
% end