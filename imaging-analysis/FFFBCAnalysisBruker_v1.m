% FFFBCAnalysisBruker_v1.m
%
% Script for analyzing binary contrast full field flash stimuli (e.g. fff or 
%  search stimulus), with option to pool

%% (1) Set up the workspace

close all;
clearvars;

% pData directory
pdataPath = '/Volumes/LEGTRACKING/L1L2_Michelle/220220/pdata_220328/';
dataPath = '/Volumes/LEGTRACKING/L1L2_Michelle/AnalyzedData';
figPath = '/Volumes/LEGTRACKING/L1L2_Michelle/Figures';

%% Read metadata spreadsheet and find indices 
% struct containing arrays or cell arrays for each column in your 
% experimental metadata spreadsheet 
readDatabase_selectSamples_NatStim;

%% (3) Select time series of interest 

% REQUIRED VARIABLES (see sections below for examples)
%   timeseriesInd
%   refStimCode (set to '' if your experiment only has one stimulus)
%   refPlotTitle
%   nrefPlotTitle (e.g. title for FFFOntoGray plots)
%   yScale - y-axis limit in your plot
%   inv  - set to 1 for ASAP voltage imaging, 0 otherwise
%   binWidthMult - used when averaging responses across trials. Determines
%       the width of the bin (in multiples of 1/frameRate) that is
%       averaged. This is a trade off between too much smoothing and too
%       much noise from limited sampling 
%   interpFrameRate - what frame rate to resample data to
%   pairedEpochs - for FFFOntoGray only to indicate light/dark epoch pairs
%   saveName - if you want the code to automatically save your data matrices 

%% L1 iGluSnFR and GECO1a test data, layer 1
timeseriesInd = find(i_L1_glusnfr_geco1a .*i_M1 .* (~metaDat.isMoving) .*...
    (i_300msSearchStimFlashRGB + i_300msFFF_RGB));

refStimCode = '300ms_searchStimFlash_RGB';
nrefStimCode = {'300msFFF_binarycontrast_RGB'};

refColumn = 1;
nrefColumn = 2;
interpFrameRate = floor(120); % sample to 120 Hz
binWidthMult = 1;

pairedEpochs = [1,2]; 
yScale = [-0.1 0.1];

% how many trials of nref stim presented for FOV
trialsPerROI = 3; 

refPlotTitle = '300ms_searchStimFlash_RGB';
nrefPlotTitle = '300msFFF_binarycontrast_RGB';
saveName = '300msFFF_binarycontrast_RGB';
inv = 0;


%% (4) Load ROI data, doing ROI-matching across time series if necessary

% If no reference stim, create an R length vector matrix where R is the 
%  number of ROIs, with info for each ROI
if isempty(refStimCode)
    roiDataMat = createROIMatrix(timeseriesInd, metaDat, pdataPath);
    
else % reference stim exists 
    % match ROIs across different time series
    roiMetaMat = matchROIsAcrossStimuli(timeseriesInd, refStimCode, metaDat);

    % for each ROI, load the pData and save the experimental metadata into 
    %  the struct. 
    roiDataMat = loadROIData(roiMetaMat, metaDat, pdataPath);
   
    % sort such that each column represents 1 stimulus
    [refColumn, roiDataMat] = sortROIDataMatByStim(roiDataMat,...
        refStimCode,nrefStimCode);
end 

%% Save data matrix (optional)
save([savePath '/' saveName], 'roiDataMat');

%% (1) Get mean stimulus-locked response for each ROI, select responding ROIs
% pooled across multiple trials
if isempty(refStimCode) % if ref doesn't exist
    [roiDataMatMeans, iResp, iInv, ~] = ...
        filterROIs_stimLockedPlots_FFFBC_noRef_pool(roiDataMat, ...
        interpFrameRate, inv, yScale, binWidthMult);
else % if reference exists
    [roiDataMatMeans iResp, iInv, ~] = ...
        filterROIs_stimLockedPlots_FFFBC_pool(roiDataMat, refColumn, ...
        interpFrameRate, inv, yScale, binWidthMult);
end 

respROIMat = roiDataMatMeans(iResp,:);

%% (2) plot responses averaged over ROIs
if size(roiDataMatMeans,2) ~= 1
    % plot reference stimulus
    fffDuration = roiDataMatMeans(1,refColumn).stimDat.Duration{1}; 
    % only responding cells
    plot_FullFieldFlash(respROIMat(:,refColumn), inv, yScale, ...
        1/interpFrameRate, fffDuration, refPlotTitle, 1);
    refFig = gcf;
    % all cells
    plot_FullFieldFlash(roiDataMatMeans(:,refColumn), inv, yScale, ...
        1/interpFrameRate, fffDuration, refPlotTitle, 1);
end 

% plot non-reference stimulus
fffDur = roiDataMatMeans(1,nrefColumn).stimDat.Duration{1};
% all cells
plot_FullFieldFlash(roiDataMatMeans(:,nrefColumn), inv, yScale, ...
    1/interpFrameRate, fffDur, nrefPlotTitle, 1);
% only responding cells
plot_FullFieldFlash(respROIMat(:,nrefColumn), inv, yScale, ...
    1/interpFrameRate, fffDur, nrefPlotTitle, 1);
nrefFig = gcf;

%% (3) Save Data
save([savePath filesep saveName], 'roiDataMat', 'roiMetaMat', ...
    'roiDataMatMeans', 'iResp', 'iInv','binWidthMult',...
    'interpFrameRate','inv','yScale', '-v7.3');
% save figures
if size(roiDataMat,2) ~= 1
    saveas(refFig,[figPath filesep saveName '_ref'],'fig');
end 
saveas(nrefFig,[figPath filesep saveName '_fffgFig'],'fig');
