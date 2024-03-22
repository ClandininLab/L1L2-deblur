% analyze_NatStim.m
%
% Script to analyze Naturalistic stimulus. Converts
%  output of selectTimeSeries.m to mean stimulus-locked response for each
%  ROI and generates some plots. Saves this analyzed data.
% Run immediately after selectTimeSeries.m without clearing the workspace
%  variables.
% Works on single time series per ROI where the Naturalistic stimulus
%  was presented or on multiple time series per ROI where the reference
%  stimulus was FullFieldFlash_BinaryContrast or SearchStimulusFlash and
%  the non-reference stimuli are all NaturalisticStimulus.
% 
% Updated: 12/11/21 - HHY
%

%% Screen responses of individual ROIs, label as responding, 
% non-responding, or inverted
if singleCol % if only one column exists
    [roiDataMat, iResp, iInv] = ...
        filterROIs_stimLockedPlots_NatStim_noRef(roiDataMat,  ...
        interpFrameRate, binWidthMult);
else % if there's more than one column
    [roiDataMat, iResp, iInv, framesPerLDCycle] = ...
        filterROIs_stimLockedPlots_NatStim(roiDataMat, refColumn, ...
        interpFrameRate, binWidthMult);
end 

respROIMat = roiDataMat(iResp,:);

%% plot responses averaged over ROIs
% plot reference stimulus
if ~singleCol
    fffDuration = roiDataMat(1,1).stimDat.Duration{1}; 
    plot_FullFieldFlash(respROIMat(:,refColumn), inv, yScale, ...
        1/interpFrameRate, fffDuration, refPlotTitle, 1);
%     plot_FullFieldFlash(roiDataMat(:,refColumn), inv, yScale, ...
%         1/interpFrameRate, fffDuration, refPlotTitle, 1);
    refFig = gcf;
end 

for i = 1:length(nrefColumn)

%     plot_NatStimFFF(roiDataMat(:,i+1), inv, yScale, nrefPlotTitle);
    plot_NatStimFFF(respROIMat(:,i), inv, yScale, nrefPlotTitle);
    nrefFig(i) = gcf;
end

%% Save Data
save([dataPath filesep saveName], 'roiMetaMat','roiDataMat', 'iResp', ...
    'iInv','binWidthMult','inv','yScale', '-v7.3');
% save figures
if ~singleCol
    saveas(refFig,[figPath filesep saveName '_ref'],'fig');
end 
for i = 1:length(nrefFig)
    saveas(nrefFig(i),[figPath filesep saveName '_natStimFig' num2str(i)],'fig');
end