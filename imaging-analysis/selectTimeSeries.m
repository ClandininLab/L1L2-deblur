% selectTimeSeries.m
%
% Script for analyzing specific experiment. Defines time series to analyze.
% Sets up for analyze_ functions for specific stimuli.

%% (1) Set up the workspace

close all;
clear all;

addMyPaths;

savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170227';
figPath = '/Users/hyang/Documents/New Imaging Analysis/Figures/170227';

%% (2) Read metadata spreadsheet and find indices 
% struct containing arrays or cell arrays for each column in your 
% experimental metadata spreadsheet 
readDatabase_selectSamples;

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


%% L2 (21Dhh) GCaMP6f
% template for single stimulus
timeseriesInd = find(i_21Dhh_gcamp6f.*i_M2.*i_2sFFF);
refStimCode = '';
refPlotTitle = 'L2 GCaMP6f 2sFFF';
saveName = 'L2_GCaMP6f_2sFFF';
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
inv = 0; 
yScale = [-1 1];

%% L2 (R53G02AD-R29G11DBD) ASAP2f para Flp Impulse Flashes (FFFoG)
% template for reference stimulus + another stimulus, FFFoG stimulus
% note the addition of nrefStimCode, nrefPlotTitle, and pairedEpochs (for
%  FFFoG stimulus)
timeseriesInd = find(i_R53G02AD_R29G11DBD_ASAP2f_para_flp...
    .*~isMoving.*isResponding.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms200ms_Gray750ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms200ms_Gray750ms';
pairedEpochs = [1,2; 3,4]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f para Flp male';
nrefPlotTitle = 'L2 ASAP2f para Flp LightDarkFlash20ms200ms Gray750ms';
saveName = 'L2_ASAP2f_para_Flp_Impulse_170306';
inv = 1;
yScale = [-0.1 0.1];

%% L2 (R53G02AD-R29G11DBD) ASAP2f para Flp White Noise
% template for white noise stimulus
timeseriesInd = find(i_R53G02AD_R29G11DBD_ASAP2f_para_flp...
    .*i_M2.*isResponding.*~isMoving.*...
    (i_300msSearchStimFlash + i_whitenoise_1D_rand9_20ms));
refStimCode = '300ms_searchStimFlash';
binWidthMult = 1;
interpFrameRate = floor(100);
refPlotTitle = '300msSearchStimFlash L2 ASAP2f para Flp male';
saveName = 'L2_ASAP2f_para_Flp_WhiteNoise_170306';
inv = 1;
yScale = [-0.1 0.1];

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
