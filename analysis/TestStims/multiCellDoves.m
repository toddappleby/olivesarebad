%% this is for loading images

cellName = '20200409Bc3_OFFSmooth';

load(['DovesData',cellName,'_Doves.mat'])

% for k=1:length(stimIndex)

fileId = fopen(['E:\Data Analysis_2020\code\Manookin Repository\manookin-package\resources\doves\images\',stimIndex{2,5}],'rb','ieee-be');
img = fread(fileId, [1536 1024], 'uint16');
fclose(fileId);
img = double(img');
img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level
backgroundIntensity = mean(img(:));%set the mean to the mean over the image
img = img.*255; %rescale s.t. brightest point is maximum monitor level
imageMatrix = uint8(img);

% subplot(2,1,1)
% set(gcf,'position',[0,0,1280,920])
% axes('Position', [.5 1 .5 .5]);
colormap(gray(256));
imagesc(imageMatrix);
set(gca,'visible','off')
% subplot(2,1,2)  
% axes('Position', [0.1 0.2 0.3 0.5]);

% plot(stimIndex{3,1})
% end

%% for data
cellName = '20200409Bc3_OFFSmooth';

load(['DovesData',cellName,'_Doves.mat'])
count = 0;
for t = 2:2:length(stimIndex)
    count = count +1;
    subplot(4,1,count)
    plot(stimIndex{3,t})
    title(t)
    ylabel("OFF Smooth",'FontWeight','bold')
end

%%
  load('E:\Data Analysis_2020\code\Manookin Repository\manookin-package\resources\dovesFEMstims20160826.mat')
  UGH = FEMdata(stimIndex{2,1}).eyeX;
  moveTraj = zeros(length(stimIndex),length(UGH));
  for t = 1:length(stimIndex)
      eyeX = FEMdata(stimIndex{1,t}).eyeX;
      eyeY = FEMdata(stimIndex{1,t}).eyeY;
     
      for z = 1:length(diff(eyeX))
      diffX = diff(eyeX);
      diffY = diff(eyeY);

      moveTraj(t,z) = moveTraj(t,z) + sqrt((diffX(z)^2) + (diffY(z)^2));
      end
  end
  %%
  
plot(moveTraj(1,:))
set(gca,'ytick',[])
title('trajectory, img ind: 2')
save('E:\Data Analysis_2020\weeklymeeting_529\trajectory.mat','moveTraj')