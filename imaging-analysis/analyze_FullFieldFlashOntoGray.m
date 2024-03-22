% analyze_FullFieldFlashOntoGray.m
%
% Script to analyze Full Field Flash Onto Gray stimulus. Converts
%  output of selectTimeSeries.m to mean stimulus-locked response for each
%  ROI and generates some plots. Saves this analyzed data.
% Run immediately after selectTimeSeries.m without clearing the workspace
%  variables.
% Works on single time series per ROI where the Full Field Flash Onto Gray
%  was presented or on multiple time series per ROI where the reference
%  stimulus was FullFieldFlash_BinaryContrast or SearchStimulusFlash and
%  the non-reference stimuli are all FullFieldFlashOntoGray.
% 
% Updated: 3/24/17 - HHY
%

%% (1) Get mean stimulus-locked response for each ROI, select responding ROIs
if size(roiDataMat,2) == 1 % if only one column exists
    [roiDataMat, iResp, iInv, framesPerLDCycle] = ...
        filterROIs_stimLockedPlots_FFFoG_noRef(roiDataMat,pairedEpochs,...
        interpFrameRate, inv, yScale, binWidthMult);
else % if there's more than one column
    [roiDataMat, iResp, iInv, framesPerLDCycle] = ...
        filterROIs_stimLockedPlots_FFFoG(roiDataMat, refColumn, ...
        pairedEpochs, interpFrameRate, inv, yScale, binWidthMult);
end 

respROIMat = roiDataMat(iResp,:);

%% (2) plot responses averaged over ROIs
if size(roiDataMat,2) ~= 1
    % plot reference stimulus
    fffDuration = roiDataMat(1,1).stimDat.Duration{1}; 
    % only responding cells
    plot_FullFieldFlash(respROIMat(:,1), inv, yScale, ...
        1/interpFrameRate, fffDuration, refPlotTitle, 1);
    refFig = gcf;
    % all cells
    plot_FullFieldFlash(roiDataMat(:,1), inv, yScale, ...
        1/interpFrameRate, fffDuration, refPlotTitle, 1);
end 

% plot non-reference stimulus
nrefColumn = 2; % which column of data to plot
% all cells
plot_FFFOntoGray(roiDataMat(:,nrefColumn), inv, yScale, nrefPlotTitle, ...
    pairedEpochs);
% only responding cells
plot_FFFOntoGray(respROIMat(:,nrefColumn), inv, yScale, nrefPlotTitle, ...
    pairedEpochs);
nrefFig = gcf;

%% (3) Save Data
save([savePath filesep saveName], 'roiDataMat', 'roiMetaMat', 'iResp', ...
    'iInv','binWidthMult','interpFrameRate','inv','yScale', '-v7.3');
% save figures
if size(roiDataMat,2) ~= 1
    saveas(refFig,[figPath filesep saveName '_ref'],'fig');
end 
saveas(nrefFig,[figPath filesep saveName '_fffgFig'],'fig');