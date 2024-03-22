% analyze_WhiteNoise_1D.m
% 
% This script analysis responses from 1D white noise stimuli 
% (full field uniform contrast, N intensity values). 
% It computes for each cell the fluorescence-weighted average (FWA)
% stimulus of length tau seconds. FWA is an estimate of the cell's 
% linear filter. Users should specify tau in section (2). 
% The FWA of tau seconds in the past and tau in the future is also computed
% in order detect any temporal offset in the FWA. Section (4) aligns
% the FWAs of all cells temporally to the FWA of a reference cell.  
% A linear prediction of each cell's response is computed based on its FWA.
% The linear prediction is then plotted against the cell's true response. 
%
% Run sections 1-6 one at a time. You can skip the optional sections.
% 
% For further details about white noise analysis, see Chichilnisky et al,
% 2000. 
% 
% IMPORTANT VARIABLES:
%   respROIMat - data matrix of responding ROIs from stimPlots.m
%   filters - FWA stimulus vector of length tau (seconds) for all cells
%   pfFilters - FWA of tau/2 in the past and tau/2s in the future for all
%       of the cells
% 
% Last updated 3/23/17 MX

%% (1) Find the white noise column in the data matrix

if length(respROIMat(1,:)) == 1 % only one column
    refColumn = 1;
    wnColumn = 1;
else % more than one column 
    for s = length(respROIMat(1,:))
        if strcmp(respROIMat(1,s).stimcode, refStimCode)
            refColumn = s;
        elseif isa(respROIMat(1,s).stimDat.obj, 'WhiteNoise_1D_FullFieldFlash')
            wnColumn = s;
        end 
    end 
end 
% wnROIs - column vector of white noise ROIs from the ROI data matrix 
wnROIs = respROIMat(:, wnColumn);

%% Play with photodiode thresholding (optional)

 r = 1;
    [stimvalIF, stimEpochTimesIF, rcStim, rcStimInd, stimEpochStartTimes, ...
    lightStartTimes, darkStartTimes, pdThresh] ...
    = processPD(wnROIs(r).stimDat, wnROIs(r).imFrameStartTimes);

    wnROIs(r).pStimDat.stimvalIF = stimvalIF;
    wnROIs(r).pStimDat.stimEpochTimesIF = stimEpochTimesIF;
    wnROIs(r).pStimDat.rcStim = rcStim;
    wnROIs(r).pStimDat.rcStimInd = rcStimInd; 
    wnROIs(r).pStimDat.stimEpochStartTimes = stimEpochStartTimes;
    wnROIs(r).pStimDat.lightStartTimes = lightStartTimes;
    wnROIs(r).pStimDat.darkStartTimes = darkStartTimes;
    wnROIs(r).pStimDat.pdThresh = pdThresh;


%% Plot stimulus sample and dFF (optional)
figure;
for r = 1:length(wnROIs)
% for r = 1
    stimVals = wnROIs(r).pStimDat.rcStim; % 1D sequence of stim contrast values
    oResp = wnROIs(r).dFF; % observed response
    
    subplot(2,1,2);
    plot(oResp); hold on;
    xlim([1000 1200]);
    ylim([-1 1]);
    ylabel('dF/F'); xlabel('stimulus frame');
    
    subplot(2,1,1);
    plot(stimVals, 'k'); hold on;
    xlim([1000 1200]);
    ylim([-.25 1.25])
    ylabel('contrast'); xlabel('stimulus frame');
end 

%% Another way of plotting (optional)
figure;
for r = 16
stimVals = wnROIs(r).pStimDat.rcStim; 
stimTimes = wnROIs(r).pStimDat.stimEpochStartTimes;
respVals = wnROIs(r).dFF; % observed response
respTimes = wnROIs(r).imFrameStartTimes;


subplot(2,1,1);
plot(stimTimes, stimVals); hold on;
ylim([-.25 1.25])
ylabel('contrast'); xlabel('time');

subplot(2,1,2);
plot(respTimes, respVals); hold on;
ylim([-1 1]);
ylabel('dF/F'); xlabel('time');
end 
    
%% Visualize Stimulus and Imaging timestsamps (optional)
r = 9;
    stimVals = wnROIs(r).pStimDat.rcStim;
    stimTimes = wnROIs(r).pStimDat.stimEpochStartTimes;
    oRespVals = wnROIs(r).dFF;
    oRespTimes = wnROIs(r).imFrameStartTimes;
    
    pdTime = wnROIs(r).stimDat.obj.Out.pdTime;
    pdData = wnROIs(r).stimDat.obj.Out.pdData;
    
    figure;
    plot(pdTime, pdData);
    hold on;
    for i = 1:100
        line([stimTimes(i) stimTimes(i)], [0 .1], ...
            'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r'); hold on;
        line([oRespTimes(i) oRespTimes(i)], [0, .1], ...
            'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');
    end

%% (2) Linear Filter Calculation via Fluorescence-Weighted Average

% stimulus window (seconds)
tau = 0.5; 
% stimulus sampling rate, aka filter resolution - try 100-200Hz
stimSampleRate = 100; % Hz

stimSampleIFI = 1/stimSampleRate; % stimulus interframe interval
winLength = tau/stimIFI; % number of frames per stim window

% initialize array to hold linear filters for each cell 
FWA = zeros(winLength, length(wnROIs));
pfFilters = zeros(winLength*2, length(wnROIs));

tic
% loops over all ROIs and computes the FWA and past-future FWA
for r = 1:length(wnROIs)
% for r = 1
% for r=[1 3 9 11 16 18]
    % fetch stimulus and response 
    stimVals = wnROIs(r).pStimDat.rcStim; % stimulus values
    stimTimes = wnROIs(r).pStimDat.stimEpochStartTimes; % stimulus transition times
    oRespVals = wnROIs(r).dFF; % observed response values
    oRespTimes = wnROIs(r).imFrameStartTimes; % observed response times
    
    % Assume perfect stimulus timing:  start and end according to times
    % reported by photodiode but assume stimulus always gets updated at the
    % same rate:
    starttime = wnROIs(r).pStimDat.stimEpochStartTimes(1);
    endtime = wnROIs(r).pStimDat.stimEpochStartTimes(end);
    stimIFI = wnROIs(r).stimDat.obj.IFI;
    stimTimes = starttime:stimIFI:endtime; % note that rawStim samples at 
    stimVals = wnROIs(r).stimDat.rawStim;
    stimVals = stimVals(1:length(stimTimes));
    
    % try shifting stimulus values by 1 frame
%     start = 1;
%     stimVals = [stimVals(start:1000); stimVals(1005:15000); stimVals(15005:end)];
    
%     % try using only a fraction of the data
%     d = 2;
%     stimVals = stimVals(1:round(length(stimVals)/d));
%     stimTimes = stimTimes(1:round(length(stimTimes)/d));
%     oRespVals = oRespVals(1:round(length(oRespVals)/d));
%     oRespTimes = oRespTimes(1:round(length(oRespTimes)/d));
    

% % Use this to debug timing issues related to white noise plot:
%     stimVals(1:10)
%     stimTimes(1)
% 
% stimUpdateTime = mean(diff(stimTimes)); % stim contrast duration
% upSampleFactor = round(stimUpdateTime/stimIFI); % number of repetitions per contrast value
% stimVals = repelem(stimVals, upSampleFactor); 
% stimTimes = interp1(stimTimes, 1:(1/upSampleFactor):length(stimTimes));

% Plot the timestamp of the front of the filter given a specific dF/F value
% w = 3000;
% iMiddle = find(stimTimes < oRespTimes(w), 1, 'last'); % same indexing as in computeFWA()
% 
% % Compare timestamps for a single response value 
% figure;
%     plot(wnROIs(r).stimDat.obj.Out.pdTime, wnROIs(r).stimDat.obj.Out.pdData);
%     hold on; 
%         line([stimTimes(iMiddle) stimTimes(iMiddle)], [0 .1], ...
%             'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r'); hold on;
%         line([oRespTimes(w) oRespTimes(w)], [0, .1], ...
%             'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k'); 
%      
%     figure;
%     plot(wnROIs(r).stimDat.obj.Out.pdTime, wnROIs(r).stimDat.obj.Out.pdData);
%     hold on; 
%     for i = 1:100
%         line([stimTimes(i) stimTimes(i)], [0 .1], ...
%             'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r'); hold on;
%         line([oRespTimes(i) oRespTimes(i)], [0, .1], ...
%             'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k'); 
%     end 
%     
    
    % upsample stimulus timestamps to a higher resolution 
    % repmat
     
    % Split data into training and testing sets 
    
%     % compute linear filter
%     [filters(:,r), ~] = computeFWA(tau, stimSampleRate, stimVals, stimTimes,...
%         oRespVals, oRespTimes);
    
    % Compute past-future filters 
    [pfFilters(:,r), ~] = computeFWA_pastFuture(tau, stimSampleRate, stimVals, stimTimes,...
    oRespVals, oRespTimes, 0);
end 
toc

%% Plot single FWAs to remove bad ROIs
t = -tau+stimIFI:stimIFI:0; 
respInd = [];

figure;
for r = 1:length(wnROIs)
    plot(t, FWA(:,r), 'b', 'LineWidth', 1);
    xlabel('time prior to response (sec)'); 
    ylabel('normalized filter amplitude');
    title(sprintf('ROI %d', r));
    userInput = input('Responding? [Y/N]: ','s');
    if (strcmpi(userInput, 'y')) % yes, responding
        respInd = [respInd r];
    end  
end 

respWNROIs = wnROIs(respInd);
respFilters = FWA(:, respInd);

%% Plot the average FWA with single FWAs from all good ROIs
% filter timing in seconds 
% t = 0:stimIFI:tau-stimIFI;
t = -tau+stimIFI:stimIFI:0;

figure;
for r = 1:length(respWNROIs)
    plot(t, respFilters(:,r), 'c', 'LineWidth', 0.5);
    hold on;
end 

avgFilter = mean(respFilters, 2); 
plot(t, avgFilter, 'b', 'LineWidth', 3);
xlabel('time prior to response (sec)'); 
ylabel('normalized filter amplitude');
title(sprintf('%.2fs FWA, sampling rates: stim = %dHz, dF/F = %dHz, N cells = %d', ...
    tau, stimSampleRate, 1/wnROIs(1).imIFI, length(wnROIs)));

%% (3) Plot average past/future FWA on top of singles from all ROIs
% filter timing in seconds 
t = -tau:stimIFI:tau-stimIFI; 
    % t = 0:stimIFI:tau-stimIFI;

figure;
for r = 1:length(wnROIs)
% for r = [1 7 14 20]
    plot(t, pfFilters(:,r), 'LineWidth', 0.5);
    hold on;
end 

avgFilter = mean(pfFilters, 2); % mean of all past-future FWAs
flies = length(unique([wnROIs.flyID]'));

plot(t, avgFilter, 'k', 'LineWidth', 3); 
line([0 0], ylim, ...
'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');
xlabel('time prior to response (sec)'); 
ylabel('normalized filter amplitude');
title(sprintf('%.2fs FWA, sampling rates: stim = %dHz, img = %dHz, N cells = %d, flies = %d', ...
    tau, stimSampleRate, 1/wnROIs(1).imIFI, length(wnROIs), flies));

%% (4) Align temporally offset signals by cross-correlation 
% Goal of this section is to aid visualization of filter shape when filters
% are time-shifted with respect to each other due to timing issues. 

% change these parameters as needed
startCut = 0.1; % seconds to cut from the past end of the filter
endCut = 0.4; % seconds to cut from the future end of the filter 
% index of the filter you want to use as the reference - all other filters
% will be aligned in time to this one
iRef = 9; 

alignedFilters = alignFiltersXCorr(pfFilters, iRef, tau, stimSampleIFI, ...
    startCut, endCut, inv);

% plot title
flies = length(unique([wnROIs.flyID]'));
title(sprintf('%.2fs FWA, sampling rates: stim = %dHz, img = %dHz, N cells = %d, flies = %d', ...
    tau, stimSampleRate, 1/wnROIs(1).imIFI, length(alignedFilters), flies));


%% (5) Visualize Static Nonlinearity
% WORK IN PROGRESS: BUGS STILL NEED TO BE SOLVED!
% Algorithm:
% This version uses the aligned pfFilters
% 1) Compute each cell's linear response prediction by convolving the
%       cell's stimulus with its own filter.
% 4) Plot each each cell's nonlinearity separately.  

% for r = 1:size(pfFilters_Shifted, 2)
for r = 16
    % fetch stimulus and response 
    stimVals = wnROIs(r).pStimDat.rcStim;
    stimTimes = wnROIs(r).pStimDat.stimEpochStartTimes; 
    oRespVals = wnROIs(r).dFF;
    oRespTimes = wnROIs(r).imFrameStartTimes;
    
% % Filter out noise in the response
%        figure;
%        subplot(1,2,1);
%        plot(oRespVals);
%             
%             % Filter out noise from the actual response (added 170306)
%             % butterworth filter
%             d = designfilt('bandstopiir','FilterOrder',2, ...
%            'HalfPowerFrequency1',0.03,'HalfPowerFrequency2',1, ...
%            'DesignMethod','butter','SampleRate', 1/imIFI);
%             oRespVals = filtfilt(d, oRespVals);  
% 
%         subplot(1,2,2);
%         plot(oRespVals);

    % Assume perfect stimulus timing:  start and end according to times
    % reported by photodiode but assume stimulus always gets updated at the
    % same rate:
    starttime = wnROIs(r).pStimDat.stimEpochStartTimes(1);
    endtime = wnROIs(r).pStimDat.stimEpochStartTimes(end);
    stimIFI = wnROIs(r).stimDat.obj.IFI;
    stimTimes = starttime:stimIFI:endtime; % note that rawStim samples at 
    stimVals = wnROIs(r).stimDat.rawStim;
    stimVals = stimVals(1:length(stimTimes));

    % get the upsampled stimulus (same code as in computeFWA.m)
    stimSampleIFI = 1/stimSampleRate;
    stimUpdateTime = mean(diff(stimTimes)); % stim contrast duration
    upSampleFactor = round(stimUpdateTime/stimSampleIFI); % number of repetitions per contrast value
    upStimVals = repelem(stimVals, upSampleFactor); % upsampled stimulus 
    upStimTimes = interp1(stimTimes, 1:(1/upSampleFactor):length(stimTimes));
    testStim = upStimVals;
    
    % upsample the true response vector to match sampling rate of predicted response
        %     testRespTimes = resample(oRespTimes, stimSampleRate, 1/wnROIs(1).imIFI);
    [testResp, upRespTimes] = resample(oRespVals, oRespTimes, stimSampleRate);
    
    % plot original response with upsampled response
    figure;
    plot(oRespTimes, oRespVals, 'o-', upRespTimes, testResp, '+-');
    legend('original response', 'upsampled response');
    
    % fetch this cell's linear filter
    cellFilter = alignedFilters(:,r);

    % compute this cell's linear prediction of its test response 
    lp = conv(testStim, flipud(cellFilter), 'same');
    
    % normalize each predicted linear response before averaging
    lpnorm = lp./max(lp); 
    testRespNorm = testResp./max(testResp); % normalize the actual response of this cell
    % or normalize it this way:
    lpnorm = (lp - min(lp))/(max(lp)-min(lp));
    testRespNorm = (testResp - min(testResp))/(max(testResp)-min(testResp));% try this
    
    % plot linear prediction on top of true response
    figure;
    plot(upStimTimes, lpnorm); hold on; plot(upRespTimes, testRespNorm);
    legend('linear response prediction', 'observed response');
     
    % plot linear prediction of the cell against its true response
    figure;
    plot(lpnorm(1:length(testRespNorm)), testRespNorm, 'b.');
   
  
end 


%% (6) Save Data
save([savePath '\' saveName], 'wnROIs', 'pfFilters', 'alignedFilters', ...
    'respROIMat', 'roiDataMat', 'refSLResp', 'iRefResp', 'iRefInv', ...
    'refStimCode', 'inv', 'yScale', 'binWidthMult', 'flashDuration', ...
   'framesPerCycle', 'imIFI', '-v7.3');

