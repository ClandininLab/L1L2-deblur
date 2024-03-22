% model_fft_test.m
% 
% script to test model filters for fft and quantification analysis

cd('/Users/hyang/Documents/New Imaging Analysis/Marjorie/AnalyzedData/whiteNoise/');
load('fff300ms_fffrand9_16ms_L2ASAP2f_refFiltered_filtered_analyzed_170915.mat');

interpT = 0:99;
interpT = interpT * 1/120;
interpAvgFilter = interp1(t,avgFilter,interpT);

testFilter = interpAvgFilter * -1;

lambda = 10;
initInd = 8;
initVal = testFilter(initInd);


% [initVal, initInd] = max(testFilter);

for i = initInd:length(testFilter)
    ind = i - initInd;
    testFilter(i) = initVal * exp(-1/lambda * ind);
end

p1InitInd = 2;
[peakVal, peakInd] = min(testFilter);
slope = -0.25069493568632;

for j = p1InitInd:peakInd
    ind = j - p1InitInd;
    testFilter(j) = slope * ind;
end

slope2 = 0.30096015222361;
yInt = -0.75208480705895;

for k = peakInd:initInd
    ind = k-peakInd;
    testFilter(k) = ind * slope2 + yInt;
end

figure; plot(interpT,testFilter);

sampRate = 120;     % sampling rate in Hz
n = 1024;
y = fft(testFilter,n);                    % Compute DFT of testFilter
m = abs(y);                               % Magnitude
p = unwrap(angle(y));                     % Phase

f = (0:length(y)-1)*sampRate/length(y);        % Frequency vector

p1 = 1:(floor(length(y)/2)+1);

figure;
subplot(2,1,1)
% loglog(f(p1),m(p1));
plot(f(p1),m(p1));
hold on;
% loglog(f2(p1),filtM(p1),'r');
title('Magnitude')



subplot(2,1,2)
plot(f(p1),p(p1)*180/pi)
hold on;
% loglog(f2(p1),filtP(p1)*180/pi,'r');
title('Phase')