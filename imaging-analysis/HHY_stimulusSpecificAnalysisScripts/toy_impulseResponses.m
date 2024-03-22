% toy_impulseResponses.m
%
% Script to make toy impulse responses. For L1/L2 paper, Figure 1.
%
% Updated: 11/25/17
%

clearvars
close all
clc

%% which Filter
whichFilter = 4;
titleStr = 'Not Biphasic, Change Phase 1';

%% Make model filter
ifi = 1/120;
t = 0:99;
t = t * ifi;

filter = zeros(size(t));

startInd = 2;
peak1Ind = 5; % 7, 5
peak2Ind = 7; % 15, 8

peak1Val = -0.75; %-0.75 % -0.8
peak2Val = 0.22; %0.22 %0.1

lambda = 10; % 10 

% decay of second phase
for i = peak2Ind:length(filter)
    ind = i - peak2Ind;
    filter(i) = peak2Val * exp(-1/lambda * ind);
end

% 1st phase to peak
slope1 = peak1Val / (peak1Ind - startInd);
for j = startInd:peak1Ind
    ind = j-startInd;
    filter(j) = slope1 * ind;
end

% 1st phase peak to 2nd phase peak
slope2 = (peak2Val-peak1Val)/(peak2Ind-peak1Ind);
for k = peak1Ind:peak2Ind
    ind = k - peak1Ind;
    filter(k) = slope2 * ind + peak1Val;
end

% plot filter
figure; plot(t,filter);
title(titleStr);
xlabel('Time (sec)');
ylabel('Amplitude');
xlim([0 0.520]);
ylim([-0.9 0.9]);
line([0 0.520],[0 0],'color',[0 0 0]); % x-axis line

%% normalize filter - to energy (value squared)

energy = 0;
for l = 1:length(filter)
    energy = energy + filter(l)^2;
end

normFilter = filter / energy;

% figure; plot(t, normFilter);

%% fft

sampRate = 1/ifi;     % sampling rate in Hz
n = 1024;
y = fft(normFilter,n);                    % Compute DFT of normFilter
m = abs(y);                               % Magnitude

f = (0:length(y)-1)*sampRate/length(y);        % Frequency vector

p1 = 1:(floor(length(y)/2)+1);

figure;
plot(f(p1),m(p1));
hold on;
title(titleStr)
xlim([0 40]);
ylim([0 2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');


%% compute metrics

% end of 2nd phase 0.25 sec from stim start
tEndPhase2 = 0.25;
% frame corresponding to end of second phase
endPhase2 = floor(tEndPhase2/ifi);

% framePeaks of 1st and 2nd phases
[framePeak1,framePeak2]  = computeFramePeaks(filter,1);

% first zero crossing
frameZero1 = computeFrameZero1(filter,framePeak1,1);

% second zero crossing
frameZero2 = computeFrameZero2(filter,framePeak1,1);

% tPeaks - always relative to first zero crossing
tPeak1 = t(framePeak1 - frameZero1);
tPeak2 = t(framePeak2 - frameZero1);

% area under curve of first phase in units %dF/F * sec
%  between first and second zero crossings
area1 = trapz(filter(frameZero1:frameZero2)) * 100 * ifi;


% area under curve of second phase in units %dF/F * sec
%  between second zero crossing and endPhase2 (0.25sec from stim
%  start)
% dark
area2 = trapz(filter(frameZero2:endPhase2)) * 100 * ifi;        

% ratio of 2 area of second phase/area of first phase
areaRatio = area2 / area1;

%% save metrics

filterAll(whichFilter,:) = filter;
normFilterAll(whichFilter,:) = normFilter;
fftAll(whichFilter,:) = y;
magAll(whichFilter,:) = m;

framePeak1All(whichFilter) = framePeak1;
framePeak2All(whichFilter) = framePeak2;
frameZero1All(whichFilter) = frameZero1;
frameZero2All(whichFilter) = frameZero2;
tPeak1All(whichFilter) = tPeak1;
tPeak2All(whichFilter) = tPeak2;
area1All(whichFilter) = area1;
area2All(whichFilter) = area2;
areaRatioAll(whichFilter) = areaRatio;

%% save all
savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/ToyImpulseResponses';
saveName = 'toyImpulseResponses171125.mat';

save([savePath filesep saveName], 'filterAll', 'normFilterAll','fftAll',...
    'framePeak1All', 'framePeak2All', 'frameZero1All', 'frameZero2All',...
    'tPeak1All','tPeak2All','area1All','area2All','areaRatioAll','t',...
    'ifi','f','p1','n','magAll','-v7.3');