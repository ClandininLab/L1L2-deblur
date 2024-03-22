%% selectTimeSeries.m

%% (1) Set up the workspace

close all;
clear all;

addMyPaths;

%% Users: change these path names

% % analysis repository (2p-imaging-analysis)
%     analysisPath = 'C:\Users\Marjo\Documents\GitHub\2p-imaging-analysis';
% % stimulus repository (2p-stim-code)
%     stimPath = 'C:\Users\Marjo\Documents\GitHub\2p-stim-code\stimCodeDevelopment';
% % folder containing your metadata spreadsheet 
%     metadataPath = 'D:\Marjorie\ClandininLabStanford\Imaging\ImagingRawData';
% % folder containing processed data (*_pdata.mat files)
%     pdataPath = 'D:\Marjorie\ClandininLabStanford\Imaging\processedData';
% % folders for saving analyzed data and figures (up to you to decide how to
% % organize them)
%     savePath = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData'; 
%     figPath = 'D:\Marjorie\ClandininLabStanford\Imaging\Figures'; 
% 
% addpath(analysisPath);
% addpath(stimPath);
% addpath(metadataPath);
% addpath(pdataPath);

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
%   pairedEpochs - for FFFOntoGray only to indicate light/dark epoch pairs
%   saveName - if you want the code to automatically save your data matrices 


%% L2 (21Dhh) GCaMP6f (template for users)
timeseriesInd = find(i_21Dhh_gcamp6f.*i_M2.*i_2sFFF);
refStimCode = '';
refPlotTitle = 'L2 GCaMP6f 2sFFF';
saveName = 'L2_GCaMP6f_2sFFF';
inv = 0; 
yScale = [-1 1];
binWidthMult = 1;

%% Test search stim
timeseriesInd = find(i_flyID.*~isMoving.*(i_300msSearchStimFlash));
refStimCode = '';
refPlotTitle = '300msFFF L2 ASAP2f';
binWidthMult = 1;
saveName = '300msFFF_L2_ASAP2f';
inv = 1; 
yScale = [-0.1 0.1];

%% L2 GCaMP6f Male White Noise
timeseriesInd = find(i_21Dhh_gcamp6f.*i_M2.*i_male.*(i_2s_searchStimFlash + ...
    i_whitenoise_1D_rand9_50ms));

refStimCode = '2s_searchStimFlash';
refPlotTitle = 'L2 GCaMP6f male 2s Flash Search Stim';
binWidthMult = 1;
saveName = 'L2_GCaMP6f_male_WN_rand9_50ms';
yScale = [-1 1];
inv = 0;

%% L2 (R53G02AD-R29G11DBD) ASAP2f No-Flp Control FFFoG
timeseriesInd = find(i_R53G02AD_R29G11DBD_ASAP2f_para_no_flp_control...
    .*~isMoving.*isResponding.*i_M2.*(i_300msSearchStimFlash + i_fullfield_LDflash20ms200ms_Gray750ms));

refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms200ms_Gray750ms';
pairedEpochs = [1,2; 3,4]; % indices in same row are flash types to be plotted together
binWidthMult = 1;

refPlotTitle = '300msSearchStimFlash No Flp Control';
nrefPlotTitle = 'LightDarkFlash20ms200ms Gray750ms';
saveName = 'L2_ASAP2f_para_NoFlpControl_Impulse_170306';
inv = 1;
yScale = [-0.1 0.1];

%% L2 (R53G02AD-R29G11DBD) ASAP2f No-Flp Control White Noise
timeseriesInd = find(i_R53G02AD_R29G11DBD_ASAP2f_para_no_flp_control...
    .*i_M2.*isResponding.*~isMoving.*(i_300msSearchStimFlash + i_whitenoise_1D_rand9_20ms));
refStimCode = '300ms_searchStimFlash';
binWidthMult = 1;
refPlotTitle = '300msSearchStimFlash L2 ASAP2f No-Flp Cntrl';
saveName = 'L2_ASAP2f_para_NoFlpControl_WhiteNoise_170306';
inv = 1;
yScale = [-0.1 0.1];

%% L2 (R53G02AD-R29G11DBD) ASAP2f para Flp Impulse Flashes (FFFG)
timeseriesInd = find(i_R53G02AD_R29G11DBD_ASAP2f_para_flp...
    .*~isMoving.*isResponding.*i_M2.*(i_300msSearchStimFlash + i_fullfield_LDflash20ms200ms_Gray750ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms200ms_Gray750ms';
pairedEpochs = [1,2; 3,4]; 
binWidthMult = 1;
refPlotTitle = '300msSearchStimFlash L2 ASAP2f para Flp male';
nrefPlotTitle = 'L2 ASAP2f para Flp LightDarkFlash20ms200ms Gray750ms';
saveName = 'L2_ASAP2f_para_Flp_Impulse_170306';
inv = 1;
yScale = [-0.1 0.1];

%% L2 (R53G02AD-R29G11DBD) ASAP2f para Flp White Noise
timeseriesInd = find(i_R53G02AD_R29G11DBD_ASAP2f_para_flp...
    .*i_M2.*isResponding.*~isMoving.*(i_300msSearchStimFlash + i_whitenoise_1D_rand9_20ms));
refStimCode = '300ms_searchStimFlash';
binWidthMult = 1;
refPlotTitle = '300msSearchStimFlash L2 ASAP2f para Flp male';
saveName = 'L2_ASAP2f_para_Flp_WhiteNoise_170306';
inv = 1;
yScale = [-0.1 0.1];

%% (4) Load ROI data, doing ROI-matching across time series if necessary

% If no reference stim, create an R x S matrix where R is the number ...
% of ROIs and S is the number of stimuli per ROI
if isempty(refStimCode)
    roiDataMat = createROIMatrix(timeseriesInd, metaDat, pdataPath);
    
else % reference stim exists 
    % match ROIs across different time series
    roiMetaMat = matchROIsAcrossStimuli(timeseriesInd, refStimCode, metaDat);
    
    % % view struct array field values across ROIs
    % {roiMetaMat(:,1).stimclass}';

    % for each ROI, load the pData and save the experimental metadata into the
    % struct. 
    roiDataMat = loadROIData(roiMetaMat, metaDat, pdataPath);
end 

clear roiMetaMat metaDat

%% Save data matrix (optional)
save([savePath '/' saveName], 'roiDataMat');
