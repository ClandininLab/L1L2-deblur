% NaturalisticStimAnalysisBruker_v2.m
%
% Script for looking at responses to the naturalistic stimulus 
%  (NaturalisticStimulus_1D_FullField).
% Somewhat quick and dirty, but getting there...
%
% Updated so that multiple trials presenting same stimuli to same ROI are
%  processed together 
%
% CREATED: 2/16/22 - HHY
%

clearvars
close all

% pData directory
pDataPath = '/Volumes/LEGTRACKING/L1L2_Michelle/220220/pdata_220328/';
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
nrefStimCode = '2s_1epoch-2_NatStim';

refColumn = 1;
nrefColumn = 2;
interpFrameRate = floor(300); % sample to 300 Hz
binWidthMult = 1;

% how many trials of nat stim presented for FOV
trialsPerROI = 3; 

refPlotTitle = [];
nrefPlotTitle = '2s_1epoch_NatStim';
saveName = '1epoch_NatStim';
inv = 1;
yScale = [-0.01 0.01];

%% L1 iGluSnFR and GECO1a test data, layer 5
timeseriesInd = find(i_L1_glusnfr_geco1a .*i_M5 .* (~metaDat.isMoving) .*...
    (i_300msSearchStimFlashRGB + i_2s_1epoch2_NatStim));

refStimCode = '300ms_searchStimFlash_RGB';
nrefStimCode = {'2s_1epoch-2_NatStim'};

refColumn = 1;
nrefColumn = 2;
interpFrameRate = floor(120); % sample to 120 Hz
binWidthMult = 1;

% how many trials of nat stim presented for FOV
trialsPerROI = 3; 

refPlotTitle = '300ms_searchStimFlash_RGB';
nrefPlotTitle = '2s_1epoch_NatStim';
saveName = '1epoch_NatStim';
inv = 0;
yScale = [-0.05 0.05];

%% L1 iGluSnFR and GECO1a test data, layer 1
timeseriesInd = find(i_L1_glusnfr_geco1a .*i_M1 .* (~metaDat.isMoving) .*...
    (i_300msSearchStimFlashRGB + i_2s_1epoch2_NatStim));

refStimCode = '300ms_searchStimFlash_RGB';
nrefStimCode = {'2s_1epoch-2_NatStim'};

refColumn = 1;
nrefColumn = 2;
interpFrameRate = floor(120); % sample to 120 Hz
binWidthMult = 1;

% how many trials of nat stim presented for FOV
trialsPerROI = 3; 

refPlotTitle = '300ms_searchStimFlash_RGB';
nrefPlotTitle = '2s_1epoch_NatStim';
saveName = '1epoch_NatStim';
inv = 0;
yScale = [-0.05 0.05];


%% Load ROI data, doing ROI-matching across time series if necessary

% If no reference stim, create an R x S matrix where R is the number ...
% of ROIs and S is the number of stimuli per ROI

% if only 1 stimulus per ROI
if (trialsPerROI == 1)
    roiDataMat = createROIMatrix(timeseriesInd, metaDat, pDataPath);
    
else % multiple trials per ROI, match across stimuli
    % match ROIs across different time series, ref stimulus exists
    if ~isempty(refStimCode)
        roiMetaMat = matchROIsAcrossStimuli(timeseriesInd, refStimCode, ...
            metaDat);
    else
        % no unique reference stimulus, pick a stimulus to use as
        %  reference, will use first trial for that as reference
        roiMetaMat = matchROIsAcrossStimuli_noUniRef(timeseriesInd, ...
            nrefStimCode, metaDat);
    end
    
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
% 1/4/22 - updated sortROIDataMatByStim to handle multiple trials of same
%  stimulus per cell
if size(roiDataMat,2) ~= 1
    singleCol = false;
    [refColumn, roiDataMat] = sortROIDataMatByStim(roiDataMat,...
        refStimCode,nrefStimCode);
else
    singleCol = true;
end

%% Screen responses of individual ROIs, label as responding, 
% non-responding, or inverted
if isempty(refStimCode) % if no reference stimulus
    [roiDataMatMeans, iResp, iInv] = ...
        filterROIs_stimLockedPlots_NatStim_noRef_pool(roiDataMat,  ...
        interpFrameRate, binWidthMult);
else % if there's a reference stimulus
    [roiDataMatMeans, iResp, iInv, framesPerLDCycle] = ...
        filterROIs_stimLockedPlots_NatStim_pool(roiDataMat, refColumn, ...
        interpFrameRate, binWidthMult);
end 

respROIMat = roiDataMatMeans(iResp,:);

%% plot responses averaged over ROIs
% plot reference stimulus
if ~singleCol
    fffDuration = roiDataMatMeans(1,refColumn).stimDat.Duration{1}; 
    plot_FullFieldFlash(respROIMat(:,refColumn), inv, yScale, ...
        1/interpFrameRate, fffDuration, refPlotTitle, 1);
%     plot_FullFieldFlash(roiDataMat(:,refColumn), inv, yScale, ...
%         1/interpFrameRate, fffDuration, refPlotTitle, 1);
    refFig = gcf;
end 

% for i = 2:length(nrefColumn)
% 
%     plot_NatStimFFF(roiDataMat(:,i+1), inv, yScale, nrefPlotTitle);
% %     plot_NatStimFFF(respROIMat(:,i), inv, yScale, nrefPlotTitle);
%     nrefFig(i) = gcf;
% end

    plot_NatStimFFF(respROIMat(:,nrefColumn), inv, yScale, nrefPlotTitle);

%% Save Data
save([savePath filesep saveName], 'roiDataMat', 'roiMetaMat', ...
    'roiDataMatMeans', 'iResp', 'iInv','binWidthMult',...
    'interpFrameRate','inv','yScale', '-v7.3');
% save figures
if ~singleCol
    saveas(refFig,[figPath filesep saveName '_ref'],'fig');
end 
for i = 1:length(nrefFig)
    saveas(nrefFig(i),[figPath filesep saveName '_natStimFig' num2str(i)],'fig');
end
