% L1L2paper_impulseResponseModeling_v1.m
%
% Generate model impulse responses for figures 1 and 4 of L1/L2 paper
%
% Updated: 2/8/18
%

clearvars
close all
clc

%% Filters
ifi = 1/120;
t = 0:99;
t = t * ifi;

% colormap
colormap lines
cmap = colormap;
close

% derivative filter (negatively signed?)
dFilter = zeros(size(t));
dFilter(1) = -1;
dFilter(2) = 1;

% biologically plausible derivative filter - balanced biphasic
bbFilter.startInd = 2;
bbFilter.peak1Val = -0.75;
bbFilter.filterLength = 24;
bbFilter = generate_balancedBiphasicFilter(bbFilter.startInd,...
    bbFilter.peak1Val,bbFilter.filterLength,length(t));

% biologically plausible monophasic filter (same shape as 1st phase of
% balanced biphasic)
mFilter.startInd = 2;
mFilter.peak1Val = -0.75;
mFilter.filterLength = 12;
mFilter = generate_balancedBiphasicFilter(mFilter.startInd,...
    mFilter.peak1Val,mFilter.filterLength,length(t));

% family of peaky biphasic filters that resemble L1/L2 impulse responses
pFilters(1).titleStr = 'Biphasic AreaRatio~5';
pFilters(1).color = cmap(4,:);
pFilters(1).startInd = 2;
pFilters(1).peak1Ind = 5;
pFilters(1).peak2Ind = 8; % 9
pFilters(1).peak1Val = -0.75;
pFilters(1).peak2Val = 0.70;
pFilters(1).lambda = 13;
pFilters(1).filter = generate_peakyBiphasicFilter(pFilters(1).startInd,...
    pFilters(1).peak1Ind, pFilters(1).peak2Ind, pFilters(1).peak1Val, ...
    pFilters(1).peak2Val, pFilters(1).lambda, length(t));

pFilters(2).titleStr = 'Biphasic AreaRatio~2';
pFilters(2).color = cmap(1,:);
pFilters(2).startInd = 2;
pFilters(2).peak1Ind = 5;
pFilters(2).peak2Ind = 8;
pFilters(2).peak1Val = -0.75;
pFilters(2).peak2Val = 0.35;
pFilters(2).lambda = 13;
pFilters(2).filter = generate_peakyBiphasicFilter(pFilters(2).startInd,...
    pFilters(2).peak1Ind, pFilters(2).peak2Ind, pFilters(2).peak1Val, ...
    pFilters(2).peak2Val, pFilters(2).lambda, length(t));

% pFilters(2).titleStr = 'Biphasic AreaRatio~1.5';
% pFilters(2).color = cmap(1,:);
% pFilters(2).startInd = 2;
% pFilters(2).peak1Ind = 5;
% pFilters(2).peak2Ind = 8;
% pFilters(2).peak1Val = -0.75;
% pFilters(2).peak2Val = 0.35;
% pFilters(2).lambda = 13;

pFilters(3).titleStr = 'Biphasic AreaRatio~1';
pFilters(3).color = cmap(6,:);
pFilters(3).startInd = 2;
pFilters(3).peak1Ind = 5;
pFilters(3).peak2Ind = 8;
pFilters(3).peak1Val = -0.75;
pFilters(3).peak2Val = 0.185; %0.2
pFilters(3).lambda = 13;
pFilters(3).filter = generate_peakyBiphasicFilter(pFilters(3).startInd,...
    pFilters(3).peak1Ind, pFilters(3).peak2Ind, pFilters(3).peak1Val, ...
    pFilters(3).peak2Val, pFilters(3).lambda, length(t));

pFilters(4).titleStr = 'Biphasic AreaRatio~0.4';
pFilters(4).color = cmap(5,:);
pFilters(4).startInd = 2;
pFilters(4).peak1Ind = 5;
pFilters(4).peak2Ind = 8; % 7
pFilters(4).peak1Val = -0.75;
pFilters(4).peak2Val = 0.075; %0.1
pFilters(4).lambda = 13;
pFilters(4).filter = generate_peakyBiphasicFilter(pFilters(4).startInd,...
    pFilters(4).peak1Ind, pFilters(4).peak2Ind, pFilters(4).peak1Val, ...
    pFilters(4).peak2Val, pFilters(4).lambda, length(t));

pFilters(5).titleStr = 'Monophasic';
pFilters(5).color = cmap(3,:);
pFilters(5).startInd = 2;
pFilters(5).peak1Ind = 5;
pFilters(5).peak2Ind = 8; % 7
pFilters(5).peak1Val = -0.75;
pFilters(5).peak2Val = 0;
pFilters(5).lambda = 13;
pFilters(5).filter = generate_peakyBiphasicFilter(pFilters(5).startInd,...
    pFilters(5).peak1Ind, pFilters(5).peak2Ind, pFilters(5).peak1Val, ...
    pFilters(5).peak2Val, pFilters(5).lambda, length(t));

pFilters(6).titleStr = 'Monophasic Slow';
pFilters(6).color = cmap(2,:);
pFilters(6).startInd = 2;
pFilters(6).peak1Ind = 7;
pFilters(6).peak2Ind = 15;
pFilters(6).peak1Val = -0.75;
pFilters(6).peak2Val = 0;
pFilters(6).lambda = 13;
pFilters(6).filter = generate_peakyBiphasicFilter(pFilters(6).startInd,...
    pFilters(6).peak1Ind, pFilters(6).peak2Ind, pFilters(6).peak1Val, ...
    pFilters(6).peak2Val, pFilters(6).lambda, length(t));

% filters with slow, not peaky 2nd phase


%% generate filter

for n = 1:length(pFilters)
    filter = zeros(size(t));
    
    % decay of second phase
    for i = pFilters(n).peak2Ind:length(filter)
        ind = i - pFilters(n).peak2Ind;
        filter(i) = pFilters(n).peak2Val * exp(-1/pFilters(n).lambda * ind);
    end

    % 1st phase to peak
    slope1 = pFilters(n).peak1Val / ...
        (pFilters(n).peak1Ind - pFilters(n).startInd);
    for j = pFilters(n).startInd:pFilters(n).peak1Ind
        ind = j-pFilters(n).startInd;
        filter(j) = slope1 * ind;
    end

    % 1st phase peak to 2nd phase peak
    slope2 = (pFilters(n).peak2Val-pFilters(n).peak1Val)/...
        (pFilters(n).peak2Ind-pFilters(n).peak1Ind);
    for k = pFilters(n).peak1Ind:pFilters(n).peak2Ind
        ind = k - pFilters(n).peak1Ind;
        filter(k) = slope2 * ind + pFilters(n).peak1Val;
    end
    
    % save filter
    pFilters(n).filter = filter;

end

%% plot filters
figure;
hold on
colormap lines
cmap = colormap;

for n = length(pFilters):-1:1
    plot(t,pFilters(n).filter,'color',pFilters(n).color,'linewidth',1.5);   
end

legend(fliplr(legends));
xlabel('Time (sec)');
ylabel('Amplitude');
title('Toy Impulse Responses');
xlim([0 0.520]);
ylim([-0.9 0.9]);
line([0 0.520],[0 0],'color',[0 0 0]); % x-axis line

%% compute quantification metrics

for n = 1:length(pFilters)
    filter = pFilters(n).filter;
    
    % end of 2nd phase 0.25 sec from stim start
    tEndPhase2 = 0.25;
    % frame corresponding to end of second phase
    endPhase2 = floor(tEndPhase2/ifi);

    % framePeaks of 1st and 2nd phases
    [pFilters(n).framePeak1,pFilters(n).framePeak2]  = ...
        computeFramePeaks(filter,1);

    % first zero crossing
    pFilters(n).frameZero1 = computeFrameZero1(filter,...
        pFilters(n).framePeak1,1);

    % second zero crossing
    pFilters(n).frameZero2 = computeFrameZero2(filter,...
        pFilters(n).framePeak1,1);

    % tPeaks - always relative to first zero crossing
    pFilters(n).tPeak1 = t(pFilters(n).framePeak1 - pFilters(n).frameZero1);
    pFilters(n).tPeak2 = t(pFilters(n).framePeak2 - pFilters(n).frameZero1);

    % area under curve of first phase in units %dF/F * sec
    %  between first and second zero crossings
    pFilters(n).area1 = trapz(filter(...
        pFilters(n).frameZero1:pFilters(n).frameZero2)) * 100 * ifi;


    % area under curve of second phase in units %dF/F * sec
    %  between second zero crossing and endPhase2 (0.25sec from stim
    %  start)
    % dark
    pFilters(n).area2 = trapz(filter(...
        pFilters(n).frameZero2:endPhase2)) * 100 * ifi;        

    % ratio of 2 area of second phase/area of first phase
    pFilters(n).areaRatio = pFilters(n).area2 / pFilters(n).area1;
    
    % invert area1 and areaRatio 
    pFilters(n).area1 = pFilters(n).area1 * -1;
    pFilters(n).areaRatio = pFilters(n).areaRatio * -1;

end

%% plot quantification metrics

fieldList = {'area1', 'area2', 'areaRatio', 'tPeak1', 'tPeak2'};
xPos = length(pFilters):-1:1;
yScale = [0 5; 0 8; 0 5; 0 0.05; 0 0.15];

for m = 1:length(fieldList)
    figure;
    hold on;
    for n = length(pFilters):-1:1
        plot(xPos(n),pFilters(n).(fieldList{m}),'o',...
            'MarkerEdgeColor',pFilters(n).color,'MarkerFaceColor',...
            pFilters(n).color,'MarkerSize',8);
    end
    legend(fliplr(legends),'location','northeast');
    title(fieldList{m});
    ylim(yScale(m,:));
    xlim([0 40]);
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
end

%% generate repeated flashes - 8 ms dark (-1), 50 ms gray (0)

flashParams.numRepeatFlashes = 6;
flashParams.repeatFlashDur = 0.008; % in sec
flashFrames = ceil(flashParams.repeatFlashDur/ifi);
flashParams.repeatGrayDur = 0.05; % in sec
grayFrames = ceil(flashParams.repeatGrayDur/ifi);

flashParams.repeatFlashStim = [-1*ones(1,flashFrames) zeros(1,grayFrames)];
flashParams.fullRepeatFlashStim = repmat(flashParams.repeatFlashStim,1,...
    flashParams.numRepeatFlashes);

%% convolve filter with repeat flash stimulus

for n = 1:length(pFilters)
    % response
    pFilters(n).filtRepeatFlashOut = conv(flashParams.fullRepeatFlashStim,...
        pFilters(n).filter);
    
    % half-wave rectified response
    pFilters(n).rectFiltRepeatFlashOut = pFilters(n).filtRepeatFlashOut;
    pFilters(n).rectFiltRepeatFlashOut(...
        pFilters(n).rectFiltRepeatFlashOut < 0) = 0;
end


% time points of response to repeated flash stimulus
tRepeatFlash = 0:(length(pFilters(1).filtRepeatFlashOut)-1);
tRepeatFlash = tRepeatFlash * ifi;

%% plot filtered repeat flash responses

figure;
hold on
colormap lines
cmap = colormap;

for n = length(pFilters):-1:1
    plot(tRepeatFlash,pFilters(n).filtRepeatFlashOut,'color',...
        pFilters(n).color,'linewidth',1.5);   
end

% plot when flashes occured
flashInd = find(flashParams.fullRepeatFlashStim ~= 0);
flashParams.repeatFlashTimes = tRepeatFlash(flashInd);
flashParams.repeatFlashVals = ones(size(flashParams.repeatFlashTimes))*1.5;

plot(flashParams.repeatFlashTimes,flashParams.repeatFlashVals,'v',...
    'MarkerEdgeColor','k','MarkerFaceColor','k');

legend(fliplr(legends),'location','southwest');
xlabel('Time (sec)');
ylabel('Response Amplitude');
title('Modeled Responses to Repeated Flash Stimulus');
xlim([0 0.3]);
ylim([-2 1.5]);
line([0 0.5],[0 0],'color',[0 0 0]); % x-axis line

set(gcf,'Position',[400 150 350 500]);

%% plot rectified filtered repeat flash responses

figure;
hold on
colormap lines
cmap = colormap;

for n = length(pFilters):-1:1
    plot(tRepeatFlash,pFilters(n).rectFiltRepeatFlashOut,'color',...
        pFilters(n).color,'linewidth',1.5);   
end

% plot when flashes occured
flashInd = find(flashParams.fullRepeatFlashStim ~= 0);
flashParams.repeatFlashTimes = tRepeatFlash(flashInd);
flashParams.repeatFlashVals = ones(size(flashParams.repeatFlashTimes))*1.5;

plot(flashParams.repeatFlashTimes,flashParams.repeatFlashVals,'v',...
    'MarkerEdgeColor','k','MarkerFaceColor','k');

legend(fliplr(legends),'location','southwest');
xlabel('Time (sec)');
ylabel('Response Amplitude');
title('Modeled Responses to Repeated Flash Stimulus, Rectified');
xlim([0 0.3]);
ylim([-2 1.5]);
line([0 0.5],[0 0],'color',[0 0 0]); % x-axis line

set(gcf,'Position',[400 150 350 500]);

%% generate alternating flashes  
% 8 ms dark (-1), 50 ms gray (0), 8 ms light (1), 50 ms gray

flashParams.numAltFlashes = 3;
flashParams.altFlashDur = 0.008; % in sec
flashFrames = ceil(flashParams.altFlashDur/ifi);
flashParams.altGrayDur = 0.05; % in sec
grayFrames = ceil(flashParams.altGrayDur/ifi);

flashParams.altFlashStim = [-1*ones(1,flashFrames) zeros(1,grayFrames)...
    ones(1,flashFrames) zeros(1,grayFrames)];
flashParams.fullAltFlashStim = repmat(flashParams.altFlashStim,1,...
    flashParams.numAltFlashes);

%% convolve filter with alternating flash stimulus

for n = 1:length(pFilters)
    % response
    pFilters(n).filtAltFlashOut = conv(flashParams.fullAltFlashStim,...
        pFilters(n).filter);
    
    % half-wave rectified response
    pFilters(n).rectFiltAltFlashOut = pFilters(n).filtAltFlashOut;
    pFilters(n).rectFiltAltFlashOut(...
        pFilters(n).rectFiltAltFlashOut < 0) = 0;
end

% time points of response to alternating flash stimulus
tAltFlash = 0:(length(pFilters(1).filtAltFlashOut)-1);
tAltFlash = tAltFlash * ifi;

%% plot filtered alternating flash responses

figure;
hold on
colormap lines
cmap = colormap;

for n = length(pFilters):-1:1
    plot(tAltFlash,pFilters(n).filtAltFlashOut,'color',pFilters(n).color,...
        'linewidth',1.5);   
end

% plot when light flashes occured
flashInd = find(flashParams.fullAltFlashStim > 0);
flashParams.altLightFlashTimes = tAltFlash(flashInd);
flashParams.altLightFlashVals = ...
    ones(size(flashParams.altLightFlashTimes))*1.5;
plot(flashParams.altLightFlashTimes,flashParams.altLightFlashVals,'vk');

% plot when dark flashes occured
flashInd = find(flashParams.fullAltFlashStim < 0);
flashParams.altDarkFlashTimes = tAltFlash(flashInd);
flashParams.altDarkFlashVals = ...
    ones(size(flashParams.altDarkFlashTimes))*1.5;
plot(flashParams.altDarkFlashTimes,flashParams.altDarkFlashVals,'v',...
    'MarkerEdgeColor','k','MarkerFaceColor','k');

legend(fliplr(legends),'location','southeast');
xlabel('Time (sec)');
ylabel('Response Amplitude');
title('Modeled Responses to Alternating Flash Stimulus');
xlim([0 0.3]);
ylim([-2 1.5]);
line([0 0.5],[0 0],'color',[0 0 0]); % x-axis line

set(gcf,'Position',[400 150 350 500]);

%% plot rectified filtered alternating flash responses

figure;
hold on
colormap lines
cmap = colormap;

for n = length(pFilters):-1:1
    plot(tAltFlash,pFilters(n).rectFiltAltFlashOut,'color',pFilters(n).color,...
        'linewidth',1.5);   
end

% plot when light flashes occured
flashInd = find(flashParams.fullAltFlashStim > 0);
flashParams.altLightFlashTimes = tAltFlash(flashInd);
flashParams.altLightFlashVals = ...
    ones(size(flashParams.altLightFlashTimes))*1.5;
plot(flashParams.altLightFlashTimes,flashParams.altLightFlashVals,'vk');

% plot when dark flashes occured
flashInd = find(flashParams.fullAltFlashStim < 0);
flashParams.altDarkFlashTimes = tAltFlash(flashInd);
flashParams.altDarkFlashVals = ...
    ones(size(flashParams.altDarkFlashTimes))*1.5;
plot(flashParams.altDarkFlashTimes,flashParams.altDarkFlashVals,'v',...
    'MarkerEdgeColor','k','MarkerFaceColor','k');

legend(fliplr(legends),'location','southeast');
xlabel('Time (sec)');
ylabel('Response Amplitude');
title('Modeled Responses to Alternating Flash Stimulus, Rectified');
xlim([0 0.3]);
ylim([-2 1.5]);
line([0 0.5],[0 0],'color',[0 0 0]); % x-axis line

set(gcf,'Position',[400 150 350 500]);

%% normalize filter - to energy (value squared)

for n = 1:length(pFilters)
    filter = pFilters(n).filter;
    energy = 0;
    for l = 1:length(filter)
        energy = energy + filter(l)^2;
    end
    
    pFilters(n).normFilter = filter / energy;
end

%% compute fft

sampRate = 1/ifi;     % sampling rate in Hz
numSamp = 1024;

for n = 1:length(pFilters)
    y = fft(pFilters(n).normFilter,numSamp);   % Compute DFT of normFilter
%     y = fft(filters(n).filter,numSamp);
    pFilters(n).fftMag = abs(y); 
    % Magnitude
    % Frequency vector
    pFilters(n).freq = (0:length(y)-1)*sampRate/length(y);   
  
    pFilters(n).p1 = 1:(floor(length(y)/2)+1);
end

%% plot frequency curves
figure;
hold on
colormap lines
cmap = colormap;

for n = length(pFilters):-1:1
    plot(pFilters(n).freq(pFilters(n).p1),pFilters(n).fftMag(pFilters(n).p1),...
        'color',pFilters(n).color,'linewidth',1.5);   
end

legend(fliplr(legends));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT of Toy Impulse Responses');
xlim([0 40]);
ylim([0 2]);

%% save all
savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/ToyImpulseResponses';
saveName = 'toyImpulseResponses180119.mat';

save([savePath filesep saveName],'filters','flashParams','t',...
    'tRepeatFlash','tAltFlash','legends','ifi','-v7.3');