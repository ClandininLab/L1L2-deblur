% analyze_FullFieldFlash_BinaryContrast.m
%
% Script to analyze Full Field Flash Binary Contrast stimulus. Converts
%  output of selectTimeSeries.m to mean stimulus-locked response for each
%  ROI and generates some plots. Saves this analyzed data.
% Run immediately after selectTimeSeries.m without clearing the workspace
%  variables.
% Works on single time series per ROI or multiple time series per ROI,
%  provided all of those time series were presented with the Full Field
%  Flash Binary Contrast or the Search Stimulus Flash stimuli
%
% Updated: 3/24/17 - HHY
%

%% (1) Get mean stimulus-locked response for each ROI, select responding ROIs

% when there are multiple stimuli presented to the same ROI (all Full Field
%  Flash Binary Contrast)
if size(roiDataMat,2) ~= 1    
    % screen responses of individual ROIs, label as responding,
    % non-responding, or inverted; compute mean response
    [roiDataMat, iRefResp, iRefInv, BIN_WIDTH] = ...
        aggregate_fff_matchedROIs(roiDataMat, refColumn, interpFrameRate, ...
        binWidthMult);
% only one stimulus per ROI
else
    flashDuration = roiDataMat(1).stimDat.Duration{1};
    
    % screen responses of individual ROIs, label as responding,
    % non-responding, or inverted; compute mean response
    [roiDataMat, iRefResp, iRefInv, framesPerCycle, BIN_WIDTH, P_VAL_THRESH] = ...
        test_aggregate(roiDataMat, flashDuration*2,interpFrameRate, 1);
end

% only the responding ROIs 
respROIMat = roiDataMat(iRefResp,:);

%% (2) Plots
% Plot average and also individual traces 
% for single stimulus per ROI
if size(roiDataMat,2) == 1  
    % plot average and individual cell responses
    plot_FullFieldFlash(respROIMat, inv, yScale, 1/interpFrameRate, ...
        flashDuration, refPlotTitle, 1);

    % Plot just the average
    plot_FullFieldFlash(refROIs, inv, yScale, 1/interpFrameRate, ...
        flashDuration, refPlotTitle, 0);
    figID_avg = gcf;
% for multiple stimuli per ROI
else
    % reference stimulus
    plot_FullFieldFlash(respROIMat(:,1),inv,yScale,1/interpFrameRate,...
        roiDataMat(1,1).stimDat.Duration{1},...
        [refPlotTitle ' ' refStimCode],0);
    figID_avg = gcf;

    % non-reference stimuli
    for i = 2:size(respROIMat,2)
        flashDuration = roiDataMat(1,i).stimDat.Duration{1};
        plot_FullFieldFlash(respROIMat(:,i), inv, yScale, ...
            1/interpFrameRate, flashDuration, ...
            [nrefPlotTitle ' ' nrefStimCode{i-1}], 0);  
    end
end
    
%% (3) Save 
% save data
if size(roiDataMat,2) == 1 
    save([savePath filesep saveName], 'respROIMat', 'roiDataMat',  ...
        'iRefResp', 'iRefInv', 'refStimCode', 'refPlotTitle','inv', ...
        'yScale', 'binWidthMult','interpFrameRate', '-v7.3');
else
    save([savePath filesep saveName], 'respROIMat', 'roiDataMat',  ...
        'iRefResp', 'iRefInv', 'refStimCode', 'refPlotTitle', ...
        'nrefStimCode', 'nrefPlotTitle', 'inv', ...
        'yScale', 'binWidthMult','interpFrameRate', '-v7.3');
end

% save figure
saveas(figID_avg,[figPath filesep saveName '_avg'],'fig');