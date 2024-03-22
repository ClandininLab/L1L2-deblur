% selectTimeSeries_L2Project_cellTypeSpecific_FFFoG.m
%
% Script for analyzing full field flash off gray, specifically for cell
%  type specific silencing experiments that are part of the L2 project.
%
% 10/21/17

%% (1) Set up the workspace
clc
close all;
clear all;

addMyPaths;

filename = '/Users/hyang/Documents/New Imaging Analysis/2p_imaging_log.txt';
savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/171021_celltypespecific';
figPath = '/Users/hyang/Documents/New Imaging Analysis/Figures/171021_celltypespecific';

%% For Marjorie's data
pDataPath = '/Users/hyang/Documents/New Imaging Analysis/Marjorie/pData';
filename = '/Users/hyang/Documents/New Imaging Analysis/Marjorie/2p_imaging_log.txt';

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


%% L2 cell type-specific silencing - UAS-TNT 
timeseriesInd = find(i_GMR16H03_asap2f_TNT...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, UAS-TNT control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, UAS-TNT control';
saveName = 'L2lexA_ASAP2f_TNT';
inv = 1;
yScale = [-0.1 0.1];

%% L2 cell type-specific silencing - R11D03AD, R19C10DBD, UAS-TNT (Lawf2)
timeseriesInd = find(i_GMR16H03_asap2f_R11D03AD_R19C10DBD_TNT...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, Lawf2 TNT experimental';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, Lawf2 experimental';
saveName = 'L2lexA_ASAP2f_Lawf2-TNT';
inv = 1;
yScale = [-0.1 0.1];

%% L2 cell type-specific silencing - R11D03AD, R19C10DBD (Lawf2)
timeseriesInd = find(i_GMR16H03_asap2f_R11D03AD_R19C10DBD...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, Lawf2-Gal4 control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, Lawf2-Gal4 control';
saveName = 'L2lexA_ASAP2f_Lawf2';
inv = 1;
yScale = [-0.1 0.1];

%% L2 cell type-specific silencing - R52H01AD, R17C11DBD, UAS-TNT (Lawf1)
timeseriesInd = find(i_GMR16H03_asap2f_R52H01AD_R17C11DBD_TNT...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, Lawf1 TNT experimental';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, Lawf1 experimental';
saveName = 'L2lexA_ASAP2f_Lawf1-TNT';
inv = 1;
yScale = [-0.1 0.1];

%% L2 cell type-specific silencing - R52H01AD, R17C11DBD (Lawf1)
timeseriesInd = find(i_GMR16H03_asap2f_R52H01AD_R17C11DBD...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, Lawf1-Gal4 control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, Lawf1-Gal4 control';
saveName = 'L2lexA_ASAP2f_Lawf1';
inv = 1;
yScale = [-0.1 0.1];

%% L2 cell type-specific silencing - R20C11AD, R48D11DBD UAS-TNT (C2/C3)
timeseriesInd = find(i_GMR16H03_asap2f_R20C11AD_R48D11DBD_TNT...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, C2/C3 TNT experimental';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, C2/C3 experimental';
saveName = 'L2lexA_ASAP2f_C2C3-TNT';
inv = 1;
yScale = [-0.1 0.1];

%% L2 cell type-specific silencing - R20C11AD, R48D11DBD (C2/C3) *
timeseriesInd = find(i_GMR16H03_asap2f_R20C11AD_R48D11DBD...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, C2/C3-Gal4 control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, C2/C3-Gal4 control';
saveName = 'L2lexA_ASAP2f_C2C3';
inv = 1;
yScale = [-0.1 0.1];

%% (4) Load ROI data, doing ROI-matching across time series if necessary
cd(pDataPath)
% If no reference stim, create an R length vector matrix where R is the 
%  number of ROIs, with info for each ROI
if isempty(refStimCode)
    roiDataMat = createROIMatrix(timeseriesInd, metaDat, pDataPath);
    
else % reference stim exists 
    % match ROIs across different time series
    roiMetaMat = matchROIsAcrossStimuli(timeseriesInd, refStimCode, metaDat);

    % for each ROI, load the pData and save the experimental metadata into 
    %  the struct. 
    roiDataMat = loadROIData(roiMetaMat, metaDat, pDataPath);
    
    % sort such that each column represents 1 stimulus
    [refColumn, roiDataMat] = sortROIDataMatByStim(roiDataMat,...
        refStimCode,nrefStimCode);
end

%% temporarily save roiDataMat, roiMetaMat
roiDataMatPart1 = roiDataMat;
roiMetaMatPart1 = roiMetaMat;

%% concatenate
roiDataMatAll = cat(1,roiDataMatPart1,roiDataMat);
roiDataMat = roiDataMatAll;

roiMetaMatAll = cat(1,roiMetaMatPart1,roiMetaMat);
roiMetaMat = roiMetaMatAll;

%% Save data matrix (optional)
save([savePath filesep saveName], 'roiDataMat', 'roiMetaMat', '-v7.3');
