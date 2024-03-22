% toy_impulseResponses_v4.m
%
% Generate toy impulse responses and compare their responses to different
%  stimuli.
%
% Updated: 2/5/18
%

clearvars
close all
clc

%% Filter parameters
ifi = 1/120;
t = 0:99;
t = t * ifi;

% colormap
colormap lines
cmap = colormap;
close

filters(1).titleStr = 'Biphasic AreaRatio~4';
filters(1).color = cmap(4,:);
filters(1).startInd = 2;
filters(1).peak1Ind = 5;
filters(1).peak2Ind = 8; % 9
filters(1).peak1Val = -0.75;
filters(1).peak2Val = 0.6;
filters(1).lambda = 13;

filters(2).titleStr = 'Biphasic AreaRatio~1.5';
filters(2).color = cmap(1,:);
filters(2).startInd = 2;
filters(2).peak1Ind = 5;
filters(2).peak2Ind = 8;
filters(2).peak1Val = -0.75;
filters(2).peak2Val = 0.3;
filters(2).lambda = 13;

filters(3).titleStr = 'Biphasic AreaRatio~0.5';
filters(3).color = cmap(5,:);
filters(3).startInd = 2;
filters(3).peak1Ind = 5;
filters(3).peak2Ind = 8; % 7
filters(3).peak1Val = -0.75;
filters(3).peak2Val = 0.095; %0.1
filters(3).lambda = 13;

filters(4).titleStr = 'Monophasic';
filters(4).color = cmap(3,:);
filters(4).startInd = 2;
filters(4).peak1Ind = 5;
filters(4).peak2Ind = 8; % 7
filters(4).peak1Val = -0.75;
filters(4).peak2Val = 0;
filters(4).lambda = 13;


% define legends
for n = 1:length(filters)
    legends{n} = filters(n).titleStr;
end

% assign filters to cell types 1: light, 2: dark
L1Filters = [4 1];
L2Filters = [3 2];

%% generate filter

for n = 1:length(filters)
    filter = zeros(size(t));
    
    % decay of second phase
    for i = filters(n).peak2Ind:length(filter)
        ind = i - filters(n).peak2Ind;
        filter(i) = filters(n).peak2Val * exp(-1/filters(n).lambda * ind);
    end

    % 1st phase to peak
    slope1 = filters(n).peak1Val / ...
        (filters(n).peak1Ind - filters(n).startInd);
    for j = filters(n).startInd:filters(n).peak1Ind
        ind = j-filters(n).startInd;
        filter(j) = slope1 * ind;
    end

    % 1st phase peak to 2nd phase peak
    slope2 = (filters(n).peak2Val-filters(n).peak1Val)/...
        (filters(n).peak2Ind-filters(n).peak1Ind);
    for k = filters(n).peak1Ind:filters(n).peak2Ind
        ind = k - filters(n).peak1Ind;
        filter(k) = slope2 * ind + filters(n).peak1Val;
    end
    
    % save filter
    filters(n).filter = filter;

end

%% plot filters
figure;
hold on
colormap lines
cmap = colormap;

for n = length(filters):-1:1
    plot(t,filters(n).filter,'color',filters(n).color,'linewidth',1.5);   
end

legend(fliplr(legends));
xlabel('Time (sec)');
ylabel('Amplitude');
title('Toy Impulse Responses');
xlim([0 0.520]);
ylim([-0.9 0.9]);
line([0 0.520],[0 0],'color',[0 0 0]); % x-axis line

%% plot filters by L1 vs L2
figure;
hold on
colormap lines
cmap = colormap;

for n = length(filters):-1:1
    plot(t,filters(n).filter,'color',filters(n).color,'linewidth',1.5);   
end

legend({'L1 light','L2 light','L2 dark','L1 dark'});
xlabel('Time (sec)');
ylabel('Amplitude');
title('Toy Impulse Responses');
xlim([0 0.520]);
ylim([-0.9 0.9]);
line([0 0.520],[0 0],'color',[0 0 0]); % x-axis line

%% compute quantification metrics

for n = 1:length(filters)
    filter = filters(n).filter;
    
    % end of 2nd phase 0.25 sec from stim start
    tEndPhase2 = 0.25;
    % frame corresponding to end of second phase
    endPhase2 = floor(tEndPhase2/ifi);

    % framePeaks of 1st and 2nd phases
    [filters(n).framePeak1,filters(n).framePeak2]  = ...
        computeFramePeaks(filter,1);

    % first zero crossing
    filters(n).frameZero1 = computeFrameZero1(filter,...
        filters(n).framePeak1,1);

    % second zero crossing
    filters(n).frameZero2 = computeFrameZero2(filter,...
        filters(n).framePeak1,1);

    % tPeaks - always relative to first zero crossing
    filters(n).tPeak1 = t(filters(n).framePeak1 - filters(n).frameZero1);
    filters(n).tPeak2 = t(filters(n).framePeak2 - filters(n).frameZero1);

    % area under curve of first phase in units %dF/F * sec
    %  between first and second zero crossings
    filters(n).area1 = trapz(filter(...
        filters(n).frameZero1:filters(n).frameZero2)) * 100 * ifi;


    % area under curve of second phase in units %dF/F * sec
    %  between second zero crossing and endPhase2 (0.25sec from stim
    %  start)
    % dark
    filters(n).area2 = trapz(filter(...
        filters(n).frameZero2:endPhase2)) * 100 * ifi;        

    % ratio of 2 area of second phase/area of first phase
    filters(n).areaRatio = filters(n).area2 / filters(n).area1;
    
    % invert area1 and areaRatio 
    filters(n).area1 = filters(n).area1 * -1;
    filters(n).areaRatio = filters(n).areaRatio * -1;

end

%% plot quantification metrics

fieldList = {'area1', 'area2', 'areaRatio', 'tPeak1', 'tPeak2'};
xPos = length(filters):-1:1;
yScale = [0 5; 0 8; 0 5; 0 0.05; 0 0.15];

for m = 1:length(fieldList)
    figure;
    hold on;
    for n = length(filters):-1:1
        plot(xPos(n),filters(n).(fieldList{m}),'o',...
            'MarkerEdgeColor',filters(n).color,'MarkerFaceColor',...
            filters(n).color,'MarkerSize',8);
    end
    legend(fliplr(legends),'location','northeast');
    title(fieldList{m});
    ylim(yScale(m,:));
    xlim([0 40]);
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
end

%% generate stimuli: pairs of 8 ms flashes, separated by 50 ms

flashParams.flashContrasts = [-1 1]; % dark and light flashes
flashParams.flashDur = 0.008; % in sec
flashFrames = ceil(flashParams.flashDur/ifi);
flashParams.grayDur = 0.05; % in sec
grayFrames = ceil(flashParams.grayDur/ifi);

counter = 1;
for i = 1:length(flashParams.flashContrasts)
    for k = 1:length(flashParams.flashContrasts)
        flashParams.flashStim(counter,:) = [...
            zeros(1,length(filters(1).filter))...
            flashParams.flashContrasts(i)*ones(1,flashFrames)...
            zeros(1,grayFrames)...
            flashParams.flashContrasts(k)*ones(1,flashFrames)...
            zeros(1,grayFrames)...
            zeros(1,length(filters(1).filter))];
        counter = counter + 1;
    end
end


%% predicted responses to stimuli

% whichFilt = L1Filters;
whichFilt = L2Filters;
% whichFilt = [1 1];

for i = 1:size(flashParams.flashStim,1)
    stim = flashParams.flashStim(i,:);
    filtLen = length(filters(whichFilt(1)).filter);
    
    resp = zeros(size(stim));
    % "convolve" with light dark weighting
    for k = 1:length(stim)
        if (stim(k) > 0) % light
            filt = filters(whichFilt(1)).filter;
            resp(k:(k+filtLen-1)) = resp(k:(k+filtLen-1)) + ...
                stim(k) * filt;
        elseif (stim(k) < 0) % dark
            filt = filters(whichFilt(2)).filter;
            resp(k:(k+filtLen-1)) = resp(k:(k+filtLen-1)) + ...
                stim(k) * filt;
        end
    end
    
    % normalize
    
    % remove zero padding
%     fullResp = resp((filtLen+1):(length(resp)-filtLen));
    fullResp = resp((filtLen+1):end);
    % save response
    flashResp(i,:) = fullResp;
end

% time points of response to flash stimulus
tFlash = 0:(length(flashResp(1,:))-1);
tFlash = tFlash * ifi;

%% plot filtered repeat flash responses

figure;
hold on
colormap lines
cmap = colormap;

for n = 1:size(flashResp,1)
    plot(tFlash,flashResp(n,:),'linewidth',1.5);   
end

% plot when flashes occured
% flashInd = find(flashParams.fullRepeatFlashStim ~= 0);
% flashParams.repeatFlashTimes = tRepeatFlash(flashInd);
% flashParams.repeatFlashVals = ones(size(flashParams.repeatFlashTimes))*1.5;
% 
% plot(flashParams.repeatFlashTimes,flashParams.repeatFlashVals,'v',...
%     'MarkerEdgeColor','k','MarkerFaceColor','k');

pairedFlashLegend = {'Dark-Dark', 'Dark-Light', 'Light-Dark', ...
    'Light-Light'};

legend(pairedFlashLegend,'location','southwest');
xlabel('Time (sec)');
ylabel('Response Amplitude');
title('Modeled L2 Responses to Paired Flash Stimulus');
xlim([0 0.3]);
ylim([-2 1.5]);
line([0 0.5],[0 0],'color',[0 0 0]); % x-axis line

set(gcf,'Position',[400 150 350 500]);

%% plot rectified filtered repeat flash responses

figure;
hold on
colormap lines
cmap = colormap;

for n = length(filters):-1:1
    plot(tRepeatFlash,filters(n).rectFiltRepeatFlashOut,'color',...
        filters(n).color,'linewidth',1.5);   
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

% flashParams.altFlashStim = [-1*ones(1,flashFrames) zeros(1,grayFrames)...
%     ones(1,flashFrames) zeros(1,grayFrames)];
flashParams.altFlashStim = [-1*ones(1,flashFrames) zeros(1,grayFrames)...
    zeros(1,flashFrames) zeros(1,grayFrames)];

flashParams.fullAltFlashStim = repmat(flashParams.altFlashStim,1,...
    flashParams.numAltFlashes);

%% convolve filter with alternating flash stimulus

for n = 1:length(filters)
    % response
    filters(n).filtAltFlashOut = conv(flashParams.fullAltFlashStim,...
        filters(n).filter);
    
    % half-wave rectified response
    filters(n).rectFiltAltFlashOut = filters(n).filtAltFlashOut;
    filters(n).rectFiltAltFlashOut(...
        filters(n).rectFiltAltFlashOut < 0) = 0;
end

% time points of response to alternating flash stimulus
tAltFlash = 0:(length(filters(1).filtAltFlashOut)-1);
tAltFlash = tAltFlash * ifi;

%% plot filtered alternating flash responses

figure;
hold on
colormap lines
cmap = colormap;

for n = length(filters):-1:1
    plot(tAltFlash,filters(n).filtAltFlashOut,'color',filters(n).color,...
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

for n = length(filters):-1:1
    plot(tAltFlash,filters(n).rectFiltAltFlashOut,'color',filters(n).color,...
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

for n = 1:length(filters)
    filter = filters(n).filter;
    energy = 0;
    for l = 1:length(filter)
        energy = energy + filter(l)^2;
    end
    
    filters(n).normFilter = filter / energy;
end

%% compute fft

sampRate = 1/ifi;     % sampling rate in Hz
numSamp = 1024;

for n = 1:length(filters)
    y = fft(filters(n).normFilter,numSamp);   % Compute DFT of normFilter
%     y = fft(filters(n).filter,numSamp);
    filters(n).fftMag = abs(y); 
    % Magnitude
    % Frequency vector
    filters(n).freq = (0:length(y)-1)*sampRate/length(y);   
  
    filters(n).p1 = 1:(floor(length(y)/2)+1);
end

%% plot frequency curves
figure;
hold on
colormap lines
cmap = colormap;

for n = length(filters):-1:1
    plot(filters(n).freq(filters(n).p1),filters(n).fftMag(filters(n).p1),...
        'color',filters(n).color,'linewidth',1.5);   
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