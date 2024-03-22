% Analysis of imaging data of FullFieldFlash stimuli 

cd('C:\Users\Marjo\Documents\GitHub\2p-stim-code\stimCodeDevelopment\timingtests\160831');

% t = Tiff('LDM_003.tif'); % change to LDM_004 or LDM_005 to view the other stimuli
% imageData = t.read;
% close(t)

% Read all frames in the time series into a Y x X x totframes matrix where
% Y and X are the dimensions of the field of view (one frame). totframes is
% the total number of frames in the time series

xdim = 200;
ydim = 20;
totframes = 3000;
I = zeros(ydim, xdim, totframes);
for i = 1:totframes
    I(:,:,i) = imread('LDM_003.tif', i);
end 

% view the images
for j = 1:100
    subplot(10,10, j);
    image(I(:,:,j+100));
end 

%% Compute average frame intensity over time. 
% Average all the pixel intensities in each frame
frAvgIntensity = zeros(totframes, 1); 

for f = 1:totframes
    framePixels = reshape(I(:,:,f), [xdim*ydim 1]);
    frAvgIntensity(f) = mean(framePixels);
end 

plot(frAvgIntensity);
xlabel('frame');
ylabel('mean pixel intensity');

%% Plot average line intensity over time

% Average line intensity over time
% lineAvgIntensity = zeros(ydim*totframes, 1);
delay = 13; % the 4.9ms frame transition is about 13 line scans 
lineAvgInt = zeros(ydim+delay, totframes);
frameTransition = zeros(delay, 1); % vector of 0s to pad the end of each frame to account for frame transition delay

for f = 1:totframes
    frameLineAvgs = mean(I(:,:,f), 2);
    lineAvgInt(:,f) = [frameLineAvgs; frameTransition];
    L = reshape(lineAvgInt, [((ydim+delay)*totframes) 1]);
end 

linescanfreq = 2800; % 1400x2 for bidirectional
time = (1:1:(ydim+delay)*totframes)./linescanfreq; % seconds

figure;
plot(time, L);
xlabel('seconds');
% xlabel('line scan');
ylabel('mean pixel intensity');
