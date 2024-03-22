% analyze_WhiteNoise_1D.m
% 
% This script analyzes responses from 1D white noise stimuli 
% (full field uniform contrast, N intensity values). For each ROI, 
% it computes the fluorescence-weighted average (FWA) stimulus of length 
% tau*2 seconds, representing the average stimulus occurring tau seconds
% prior to an observed response and tau seconds after. FWA is an estimate 
% of the cell's linear filter. Users should specify tau in Section (3).
% We compute a "past-future" FWA instead of a regular FWA up to time zero 
% because it is useful for detecting any temporal offset in the FWA, which
% may result from timing issues upstream in the analysis pipeline. 
% A prediction of the cell's linear response is computed based on its FWA.
% Next, static nonlinearity of each cell is visualized by 
% plotting the true response of the cell against its linear prediction.
% Section (4) aligns the FWAs of all cells temporally to the FWA of 
% a reference cell.  
%
% Run sections 1-6 one at a time. You can skip the optional sections.
% 
% For further details about white noise analysis, see Chichilnisky et al,
% 2000. 
% 
% IMPORTANT VARIABLES:
%   respROIMat - data matrix of responding ROIs from stimPlots.m
%   pfFilters - FWA of tau seconds in the past and tau in the future for
%       all ROIs
%   alignedFilters - matrix of FWAs aligned to a reference FWA 
%       (indexed iRef)
% 
% Last updated 3/26/17 MX

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

%% (2) Plot stimulus and dFF sample (optional)

figure;
for r = [1 2 3] % pick which ROIs you want to plot
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
    
%% (3) Linear Filter Calculation via Fluorescence-Weighted Average

% stimulus window (seconds)
tau = 0.5; 
% stimulus sampling rate, aka filter resolution - try 100-200Hz
stimSampleRate = interpFrameRate; % Hz

stimSampleIFI = 1/stimSampleRate; % stimulus interframe interval
winLength = tau/stimSampleIFI; % number of frames per stim window

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

    % Split data into training and testing sets 
    
%     % compute linear filter
%     [filters(:,r), ~] = computeFWA(tau, stimSampleRate, stimVals, stimTimes,...
%         oRespVals, oRespTimes);
    
    % Compute past-future filters 
    [pfFilters(:,r), upStimVals, upStimTimes] = computeFWA_pastFuture(tau, stimSampleRate, stimVals, stimTimes,...
    oRespVals, oRespTimes, 0);

    % ----------- Part II of white noise analysis ---------- %
    testStim = upStimVals;
    
    % upsample the true response vector to match sampling rate of predicted response
        %     testRespTimes = resample(oRespTimes, stimSampleRate, 1/wnROIs(1).imIFI);
    [testResp, upRespTimes] = resample(oRespVals, oRespTimes, stimSampleRate);
    
    % plot original response with upsampled response
    figure;
    plot(oRespTimes, oRespVals, 'o-', upRespTimes, testResp, '+-');
    legend('original response', 'upsampled response');
    
    % fetch this cell's linear filter
    cellFilter = pfFilters(:,r);

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
toc

%% (4) Plot average past/future FWA on top of singles from all ROIs
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

%% (5) Align temporally offset signals by cross-correlation (optional)
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

%% (6) Save Data
save([savePath filesep saveName], 'wnROIs', 'pfFilters', 'alignedFilters', ...
    'respROIMat', 'roiDataMat', 'refSLResp', 'iRefResp', 'iRefInv', ...
    'refStimCode', 'inv', 'yScale', 'binWidthMult', 'flashDuration', ...
   'framesPerCycle', 'imIFI', '-v7.3');

