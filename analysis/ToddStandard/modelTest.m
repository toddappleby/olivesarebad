%FILTER PARAMETERS


%starter parameters 
sdCenter= 50;
sdSurround= 75;
ampCenter=1;
ampSurround=.5;



[x y]=meshgrid(linspace(-6*50,6*50,12));
r = sqrt(x.^2+y.^2);

%stim tests

testStim = ones(12);
testStim(:,[2 4 6 8 10 12])=0;

testStim2 = ones(12);
testStim2(:,[1 3 5 7 9 11])=0;




%% 2D dog

dogRF = exp(-(r)/(2*sdCenter.^2))-ampSurround*exp(-(r)/(2*sdSurround^2));  %


%% convolve stim and rf

%match stim size w/ rf size?


% 
tester1 = testStim(:)'*dogRF(:);
% testStim2(:)'*dogRF(:);
% testStim3(:)'*dogRF(:);


% resp1 = conv2(dogRF,testStim);
% resp2 = conv2(dogRF,testStim2);
% resp3 = conv2(dogRF,testStim3);
% resp4 = conv2(dogRF,testStim4);
% resp5 = conv2(dogRF,testStimBig);

%% image loader
load('C:\Users\reals\Documents\manookin-package\resources\doves\dovesFEMstims20160826.mat')
fileId = fopen(['C:\Users\reals\Documents\manookin-package\resources\doves\images\',FEMdata(6).ImageName],'rb','ieee-be');
img = fread(fileId, [1536 1024], 'uint16');
fclose(fileId);
img = double(img');
img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level
img2 = img;
backgroundIntensity = mean(img(:));%set the mean to the mean over the image
img = img.*255; %rescale s.t. brightest point is maximum monitor level


imagesc(imageMatrix)
colormap(gray)

%% parabolas (gain)

x=-10:0.1:10;
y=sqrt((36+x.^4.1)/9);
plot(x(length(x)/2:end),y(length(y)/2:end))
hold on
y=sqrt((36+x.^4.3)/9);
plot(x(length(x)/2:end),y(length(y)/2:end))
y=sqrt((36+x.^4.5)/9);
plot(x(length(x)/2:end),y(length(y)/2:end))

%% parabolas (xshift)

x=-10:0.1:10;
y=sqrt((36+x.^4)/9);

plot(x(length(x)/2:end),y(length(y)/2:end))
hold on
x=-12:0.1:8;

plot(x(length(x)/2:end),y(length(y)/2:end))
x=-14:0.1:6;

plot(x(length(x)/2:end),y(length(y)/2:end))