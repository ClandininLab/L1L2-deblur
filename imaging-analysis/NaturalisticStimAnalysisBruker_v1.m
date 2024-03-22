% NaturalisticStimAnalysisBruker_v1.m
%
% Script for looking at responses to the naturalistic stimulus 
%  (NaturalisticStimulus_1D_FullField).
% Somewhat quick and dirty, but getting there...
%
% CREATED: 11/23/21 - HHY
%

clearvars
close all

% pData directory
pDataPath = '/Volumes/LEGTRACKING/L1L2_Michelle/pData/';
dataPath = '/Volumes/LEGTRACKING/L1L2_Michelle/AnalyzedData';
figPath = '/Volumes/LEGTRACKING/L1L2_Michelle/Figures';

%% Read metadata spreadsheet and find indices 
% struct containing arrays or cell arrays for each column in your 
% experimental metadata spreadsheet 
readDatabase_selectSamples_NatStim;

%% Select time series of interest

%% 2 epoch NatStim
% timeseriesInd = find(i_21Dhh_asap2f_jrgeco1b.*i_M2.*(~metaDat.isMoving).*...
%     (i_300msSearchStimFlash + i_2s_2epoch_NatStim).*i_920);
% 
% refStimCode = '300ms_searchStimFlash';
% nrefStimCode = '2s_2epoch_NatStim';
% 
% refColumn = 1;
% nrefColumn = 2:3;
% interpFrameRate = floor(120); % sample to 120 Hz
% binWidthMult = 1;
% 
% refPlotTitle = '300msSearchStimFlash';
% nrefPlotTitle = '2s_2epoch_NatStim';
% saveName = '2epoch_NatStim';
% inv = 1;
% yScale = [-0.08 0.08];

%% 1 epoch NatStim
% timeseriesInd = find(i_21Dhh_asap2f_jrgeco1b.*i_M2.*(~metaDat.isMoving).*...
%     (i_300msSearchStimFlash + i_2s_1epoch_NatStim).*i_920);
% 
% refStimCode = '300ms_searchStimFlash';
% nrefStimCode = '2s_1epoch_NatStim';
% 
% refColumn = 1;
% nrefColumn = 2;
% interpFrameRate = floor(120); % sample to 120 Hz
% binWidthMult = 1;
% 
% refPlotTitle = '300msSearchStimFlash';
% nrefPlotTitle = '2s_1epoch_NatStim';
% saveName = '1epoch_NatStim';
% inv = 1;
% yScale = [-0.08 0.08];

%% 1 epoch NatStim, no ref
timeseriesInd = find(i_21Dhh_asap2f.*i_M2.*(~metaDat.isMoving).*...
    (i_2s_1epoch2_NatStim).*i_920);

refStimCode = [];
nrefStimCode = '2s_1epoch_NatStim';

refColumn = 1;
nrefColumn = 2;
interpFrameRate = floor(150); % sample to 300 Hz
binWidthMult = 1;

refPlotTitle = [];
nrefPlotTitle = '2s_1epoch_NatStim';
saveName = '1epoch_NatStim';
inv = 1;
yScale = [-0.01 0.01];

%% Load ROI data, doing ROI-matching across time series if necessary

% If no reference stim, create an R x S matrix where R is the number ...
% of ROIs and S is the number of stimuli per ROI
if isempty(refStimCode)
    roiDataMat = createROIMatrix(timeseriesInd, metaDat, pDataPath);
    
else % reference stim exists 
    % match ROIs across different time series
    roiMetaMat = matchROIsAcrossStimuli(timeseriesInd, refStimCode, metaDat);
    
    % % view struct array field values across ROIs
    % {roiMetaMat(:,1).stimclass}';

    % for each ROI, load the pData and save the experimental metadata into the
    % struct. 
    roiDataMat = loadROIData(roiMetaMat, metaDat, pDataPath);
end 

%% Sort reference and non-reference
% updated: 3/16/17 - package sorting into function
% 7/15/21 - note: this doesn't work correctly on ROIs matched across
%  stimuli when the stimulus is the same across multiple trials
%  for now, just don't run this section - fix later (or restrict to 1
%  stimulus presentation per type per cell)
if size(roiDataMat,2) ~= 1
    singleCol = false;
    [refColumn, roiDataMat] = sortROIDataMatByStim(roiDataMat,...
        refStimCode,nrefStimCode);
else
    singleCol = true;
end

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
