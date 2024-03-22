% selectTimeSeries.m
%
% Script for analyzing specific experiment. Defines time series to analyze.
% Sets up for analyze_ functions for specific stimuli.

%% (1) Set up the workspace

close all;
clear all;

addMyPaths;

filename = '/Users/hyang/Documents/New Imaging Analysis/2p_imaging_log.txt';
savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170621';
figPath = '/Users/hyang/Documents/New Imaging Analysis/Figures/170621';

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


%% L2>ASAP2f (GMR16H03), 482/18, ND 1.3
timeseriesInd = find(i_GMR16H03_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms + ...
    i_fullfield_LDflash20ms_c025Gray500ms).*...
    (i_ND == 1.3).* i_482_18);
refStimCode = '300ms_searchStimFlash';
nrefStimCode = {'fullfield_LDflash20ms_Gray500ms',...
    'fullfield_LDflash20ms_c025Gray500ms'};
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f 482/18 ND 1.3';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f 482/18 ND 1.3';
saveName = 'GMR16H03_ASAP2f_482-18_ND1-3_FFFoG';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (GMR16H03), 482/18, ND 1
timeseriesInd = find(i_GMR16H03_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms + ...
    i_fullfield_LDflash20ms_c025Gray500ms).*...
    (i_ND == 1).* i_482_18);
refStimCode = '300ms_searchStimFlash';
nrefStimCode = {'fullfield_LDflash20ms_Gray500ms',...
    'fullfield_LDflash20ms_c025Gray500ms'};
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f 482/18 ND 1.0';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f 482/18 ND 1.0';
saveName = 'GMR16H03_ASAP2f_482-18_ND1_FFFoG';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (GMR16H03), 447/60, ND 2.0
timeseriesInd = find(i_GMR16H03_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms + ...
    i_fullfield_LDflash20ms_c025Gray500ms).*...
    (i_ND == 2).* i_447_60);
refStimCode = '300ms_searchStimFlash';
nrefStimCode = {'fullfield_LDflash20ms_Gray500ms',...
    'fullfield_LDflash20ms_c025Gray500ms'};
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f 447/60 ND 2.0';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f 447/60 ND 2.0';
saveName = 'GMR16H03_ASAP2f_447-60_ND2_FFFoG';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (21Dhh), 10 ms flashes
timeseriesInd = find(i_21Dhh_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash10ms_Gray500ms).*...
    (i_ND == 1.3).* i_482_18 .* (i_pwm == 200));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash10ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f';
nrefPlotTitle = 'LightDarkFlash10ms Gray500ms L2 ASAP2f';
saveName = '21Dhh_ASAP2f_FFFoG_10msFlash';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (21Dhh), 20 ms flashes
timeseriesInd = find(i_21Dhh_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms).*...
    (i_ND == 1.3).* i_482_18 .* (i_pwm == 200));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f';
saveName = '21Dhh_ASAP2f_FFFoG_20msFlash';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (21Dhh), 30 ms flashes
timeseriesInd = find(i_21Dhh_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash30ms_Gray500ms).*...
    (i_ND == 1.3).* i_482_18 .* (i_pwm == 200));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash30ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f';
nrefPlotTitle = 'LightDarkFlash30ms Gray500ms L2 ASAP2f';
saveName = '21Dhh_ASAP2f_FFFoG_30msFlash';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (21Dhh), 40 ms flashes
timeseriesInd = find(i_21Dhh_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash40ms_Gray500ms).*...
    (i_ND == 1.3).* i_482_18 .* (i_pwm == 200));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash40ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f';
nrefPlotTitle = 'LightDarkFlash40ms Gray500ms L2 ASAP2f';
saveName = '21Dhh_ASAP2f_FFFoG_40msFlash';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (21Dhh), dark (0, 0.125, 0.25 from 0.5; 0, 0.125 from 0.25)
timeseriesInd = find(i_21Dhh_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + ...
    i_fullfield_Dc0c0125c025flash20ms_c05Gray500ms + ...
    i_fullfield_Dc0c0125flash20ms_c025Gray500ms).*...
    (i_ND == 1.3).* i_482_18 .* (i_pwm == 200));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = {'fullfield_D0-0125flash20ms_025Gray500ms',...
    'fullfield_D0-0125-025flash20ms_05Gray500ms'};
pairedEpochs = [1;2;3]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f';
nrefPlotTitle = 'DarkFlashes 20ms Gray500ms L2 ASAP2f';
saveName = '21Dhh_ASAP2f_FFFoG_DarkFlash';
inv = 1;
yScale = [-0.1 0.1];

%% L2>ASAP2f (21Dhh), light (1, 0.75, 0.5 from 0.25; 1, 0.75 from 0.25)
timeseriesInd = find(i_21Dhh_asap2f.*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + ...
    i_fullfield_Lc075c1flash20ms_c05Gray500ms + ...
    i_fullfield_Lc05c075c1flash20ms_c025Gray500ms).*...
    (i_ND == 1.3).* i_482_18 .* (i_pwm == 200));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = {'fullfield_L075-1flash20ms_05Gray500ms',...
    'fullfield_L05-075-1flash20ms_025Gray500ms'};
pairedEpochs = [1;2;3]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f';
nrefPlotTitle = '20 ms Light Flashes 500 ms Gray L2 ASAP2f';
saveName = '21Dhh_ASAP2f_FFFoG_LightFlash';
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
