% selectTimeSeries.m
%
% Script for analyzing specific experiment. Defines time series to analyze.
% Sets up for analyze_ functions for specific stimuli.

%% (1) Set up the workspace

close all;
clear all;

addMyPaths;

filename = '/Users/hyang/Documents/New Imaging Analysis/2p_imaging_log.txt';
savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170601';
figPath = '/Users/hyang/Documents/New Imaging Analysis/Figures/170601';

%% For Marjorie's data
pDataPath = '/Users/hyang/Documents/New Imaging Analysis/MarjoriePData/processedData';
filename = '/Users/hyang/Documents/New Imaging Analysis/MarjoriePData/2p_imaging_log.txt';

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


%% L1>ASAP2f (GMR37E04), test
timeseriesInd = find(i_GMR37E04_asap2f.*~isMoving.*i_M1.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms + ...
    i_fullfield_LDflash20ms_c025Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = {'fullfield_LDflash20ms_Gray500ms',...
    'fullfield_LDflash20ms_c025Gray500ms'};
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L1 layer M1 ASAP2f';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L1 layer M1 ASAP2f';
saveName = 'GMR37E04_ASAP2f';
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
