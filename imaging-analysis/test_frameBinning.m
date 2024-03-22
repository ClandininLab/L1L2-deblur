% test_frameBinning.m
%
% Test image frame binning strategy on mock data
%

clear all
close all

%% Mock 300 ms fff data

sigma = 1; % for adding noise
flashTime = 300; % in ms

% Mock response to 300 ms fff, 0.001ms bins; response b/w -1 and 1, centered at 0 
x = 0:299999;
y1 = exp(log(2)*x(1:20000)/20000)-1;
y2 = exp(-(x(20001:end)-20000)/80000);
y = [y1 y2];

actResp = [y, -1*y];
x = [x, x+300000];

% x is in bins, t in ms
t = 0.001 * x;

% filter response by approximate ASAP2f kinetics
filtT = 0:0.001:10;
linFilt = 0.75*exp(-2.5*filtT)+0.25*exp(-145*filtT);
resp = conv(actResp,linFilt);
resp = resp(1:length(actResp)) / max(resp(1:length(actResp)));

figure;
hold on;
plot(t,actResp,'m');
plot(t,resp,'k');
xlabel('time (ms)');
ylabel('response');

numImgFrames = 24000;

%% 

% sample img parameters
ifi = 1/360 * 1000; % in ms
imgFrameTimes = (0:(numImgFrames-1)) * ifi;

% mock stimulus transition times
dtlTimes = (1:2:(floor(imgFrameTimes(end)/flashTime)))*flashTime;
ltdTimes = (0:2:(floor(imgFrameTimes(end)/flashTime)))*flashTime;

binFrameRate = 240; % in Hz
numFramesFlash = flashTime/1000*binFrameRate;

% mock response
imgFrameTimesRel = mod(imgFrameTimes,flashTime*2);
imgInd = zeros(1,length(imgFrameTimesRel));
for i = 1:length(imgFrameTimesRel)
    imgInd(i) = find(abs(imgFrameTimesRel(i)-t)<0.0005);
end
dFF = random('norm',resp(imgInd),sigma);

% plot scatterplot of values, relative to stimulus
figure;
scatter(imgFrameTimesRel,dFF,'.');
hold on;
plot(t,resp,'r');


% Test binning code

% for binning
darkVals = cell(numFramesFlash,1);
lightVals = darkVals;

% place each imaging frame in appropriate bin
for k = 1:length(imgFrameTimes)

    ft = imgFrameTimes(k); % current frame time

    % time difference between stimulus transitions and current
    %  frame time
    tDiffLight = ft - dtlTimes;
    tDiffDark = ft - ltdTimes;

    % remove negative values - stim transition happened after frame
    tDiffLightValid = tDiffLight(tDiffLight>=0);
    tDiffDarkValid = tDiffDark(tDiffDark>=0);

    % minimum value and index in time differences
    minLightTimeDiff = min(tDiffLightValid);
    minDarkTimeDiff= min(tDiffDarkValid);

    % put imaging frame into appropriate bin, light or dark
    % depending on which one is smaller
    % dark
    if (isempty(minLightTimeDiff) || (minDarkTimeDiff < minLightTimeDiff))
        whichBin = ceil(minDarkTimeDiff / (1/binFrameRate)/1000);
        if (whichBin == 0) % corrects for index=0
            whichBin = 1;
        end
        % when stimulus presented slightly longer than expected,
        %  ignore that img frame
        if (whichBin <= numFramesFlash)
            darkVals{whichBin} = [darkVals{whichBin} dFF(k)];
        end
    else % light
        whichBin = ceil(minLightTimeDiff / (1/binFrameRate)/1000);
        if (whichBin == 0) % corrects for index=0
            whichBin = 1;
        end
        % when stimulus presented slightly longer than expected,
        %  ignore that img frame
        if (whichBin <= numFramesFlash)
            lightVals{whichBin} = [lightVals{whichBin} dFF(k)];
        end
    end
end

meanDark = zeros(1,numFramesFlash);
meanLight = meanDark;
stdErrDark = meanDark;
stdErrLight = meanDark;

for l = 1:(numFramesFlash)
    meanDark(l) = mean(darkVals{l});
    stdErrDark(l) = std(darkVals{l})/sqrt(length(darkVals{l}));
    meanLight(l) = mean(lightVals{l});
    stdErrLight(l) = std(lightVals{l})/sqrt(length(lightVals{l}));
end

[meanDark, stdErrDark] = correctFFFBinNaN(meanDark, meanLight,...
    stdErrDark);
[meanLight, stdErrLight] = correctFFFBinNaN(meanLight, meanDark,...
    stdErrLight);

% plot binned, average responses, as if real data
% binT = (1:4*numFramesFlash)/binFrameRate*1000-(1/binFrameRate*1000);
binT = (1:4*numFramesFlash)/binFrameRate*1000-(1/binFrameRate*1000/2);

figure;
plot_err_patch_v2(binT,[meanDark meanLight meanDark meanLight],...
    [stdErrDark stdErrLight stdErrDark stdErrLight], [0 0 1],[0.5 0.5 1]);
hold on;
plot([t (t+t(end))],[resp resp],'r');


% test sliding window average binning
binWidth = 25; % in ms
binShift = 1/binFrameRate *1000; % in ms

% pre-allocate 
slidingMean = zeros(1,2*flashTime/binShift);
slidingStdErr = slidingMean;

numEdgeBins = binWidth/2/binShift; 

for j=1:(2*flashTime/binShift)
    binCenter = (j-1)*binShift;
    binStart = binCenter - binWidth/2;
    binEnd = binCenter + binWidth/2;
    
    % wrap bins around
    if (j <= numEdgeBins)
        binStart = binStart + 2*flashTime;
        valInd = find((imgFrameTimesRel<=binEnd) + ...
            (imgFrameTimesRel>binStart));
    elseif (j==(2*flashTime/binShift-numEdgeBins))
        valInd = find((imgFrameTimesRel>binStart) .* ...
            (imgFrameTimesRel<=binEnd) + (imgFrameTimesRel == 0));
    elseif (j > (2*flashTime/binShift-numEdgeBins))
        binEnd = binEnd - 2*flashTime;
        valInd = find((imgFrameTimesRel>binStart) +...
            (imgFrameTimesRel<=binEnd));
    else
        valInd = find((imgFrameTimesRel>binStart) .*...
            (imgFrameTimesRel<=binEnd));
    end
    slidingMean(j) = mean(dFF(valInd));
    slidingStdErr(j) = std(dFF(valInd))/sqrt(length(valInd));
end

% compute response binned to 120Hz
downSampResp = zeros(1,flashTime*2/1000*binFrameRate);
numFrames = floor(length(resp)/length(downSampResp));
for i=1:length(downSampResp)
    whichFrames = ((numFrames*(i-1)+1):((numFrames*i)+1));
    downSampResp(i) = mean(resp(whichFrames));
end

movAvgDownSampResp = downSampResp;
for i=1:length(downSampResp)
    if (i==1)
        movAvgDownSampResp(i) = mean([downSampResp(1) downSampResp(2)]);
    elseif (i==length(downSampResp))
        movAvgDownSampResp(i)=mean([downSampResp(i-1) downSampResp(i)]);
    else
        movAvgDownSampResp(i) = mean(downSampResp((i-1):(i+1)));
    end
end
downSampT = (0:(length(downSampResp)-1))/binFrameRate*1000+(1/binFrameRate*1000/2);

% figure;
% plot_err_patch_v2(binT,[meanDark meanLight meanDark meanLight],...
%     [stdErrDark stdErrLight stdErrDark stdErrLight], [0 0 1],[0.5 0.5 1]);
% plot_err_patch_v2((1:length(slidingMean)*2)*binShift-(binWidth/(binWidth/binShift)),...
%     [slidingMean slidingMean], [slidingStdErr slidingStdErr], [1 0 0],...
%     [1 0.5 0.5]);
% hold on;
% plot([t (t+t(end))],[resp resp],'k');

figure;
hold on;
binT = (1:2*numFramesFlash)/binFrameRate*1000-(1/binFrameRate*1000/2);
plot(t,resp,'k', 'LineWidth',2);
% plot(downSampT,downSampResp,'c','LineWidth',2);
plot(downSampT,movAvgDownSampResp,'c','LineWidth',2);
plot_err_patch_v2(binT,[meanDark meanLight],...
    [stdErrDark stdErrLight], [0 0 1],[0.5 0.5 1]);
plot_err_patch_v2((1:length(slidingMean))*binShift-(binWidth/(binWidth/binShift)),...
    slidingMean, slidingStdErr, [1 0 0],...
    [1 0.5 0.5]);
xlabel('time (ms)');
ylabel('response');
