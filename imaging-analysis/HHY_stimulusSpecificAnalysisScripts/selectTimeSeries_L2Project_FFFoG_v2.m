% selectTimeSeries_L2Project_FFFoG_v2.m
%
% Script for analyzing full field flash off gray, specifically for L2 
%  project.
%
% 10/2/17

%% (1) Set up the workspace

close all;
clear all;

addMyPaths;

filename = '/Users/hyang/Documents/New Imaging Analysis/2p_imaging_log.txt';
savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/171002';
figPath = '/Users/hyang/Documents/New Imaging Analysis/Figures/171002';

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


%% L2 ort rescue - UAS-ort, ort[1], Df(3R)BSC809 experimental
timeseriesInd = find(i_21Dhh_asap2f_UASort_ort1new_Df3RBSC809...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort rescue experimental';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort rescue experimental';
saveName = 'L2_ASAP2f_UASort_ort1_Df3RBSC809';
inv = 1;
yScale = [-0.1 0.1];

%% L2 ort rescue - UAS-ort, ort[1], control
timeseriesInd = find(i_21Dhh_asap2f_UASort_ort1new...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2];  
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f UAS-ort, ort[1] control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f UAS-ort, ort[1] control';
saveName = 'L2_ASAP2f_UASort_ort1';
inv = 1;
yScale = [-0.1 0.1];

%% L2 ort rescue - UAS-ort, Df(3R)BSC809, control
timeseriesInd = find(i_21Dhh_asap2f_UASort_Df3RBSC809...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2];  
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f UAS-ort, Df(3R)BSC809 control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f UAS-ort, Df(3R)BSC809 control';
saveName = 'L2_ASAP2f_UASort_Df3RBSC809';
inv = 1;
yScale = [-0.1 0.1];

%% L2 ort rescue - ort[1], Df(3R)BSC809 control
timeseriesInd = find(i_21Dhh_asap2f_ort1new_Df3RBSC809...
    .*~isMoving.*i_M2.*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2];  
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort[1], Df(3R)BSC809 control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort[1], Df(3R)BSC809 control';
saveName = 'L2_ASAP2f_ort1_Df3RBSC809';
inv = 1;
yScale = [-0.1 0.1];


%% L2 cell type-specific silencing, L2: GMR16H03>>ASAP2f, TNT
timeseriesInd = find(i_GMR16H03_asap2f_lexAopTNT...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, L2 TNT experimental';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, L2 TNT experimental';
saveName = 'L2lexA_ASAP2f_lexAopTNT';
inv = 1;
yScale = [-0.1 0.1];

%% L2 cell type-specific silencing, L2: GMR16H03>>ASAP2f control
timeseriesInd = find(i_GMR16H03_asap2f...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f control';
saveName = 'L2lexA_ASAP2f';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, ort[C1-3]>>TNT (21Dhh, ort-lexA), pwm 200 
timeseriesInd = find(i_21Dhh_asap2f_ort_TNT...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort TNT';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort TNT';
saveName = '21Dhh_ASAP2f_ortC1-3_TNT_200pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, ort[C1-3] control (21Dhh, ort-lexA), pwm 200 
timeseriesInd = find(i_21Dhh_asap2f_ort...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort control';
saveName = '21Dhh_ASAP2f_ortC1-3_200pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, TNT control (21Dhh, ort-lexA), pwm 200
timeseriesInd = find(i_21Dhh_asap2f_TNT...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f TNT control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f TNT control';
saveName = '21Dhh_ASAP2f_TNT_200pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, ort[C1-3]>>TNT (21Dhh, ort-lexA), pwm 20 
timeseriesInd = find(i_21Dhh_asap2f_ort_TNT...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    (i_21Dhh_asap2f_ort_TNT .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort TNT 20 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort TNT 20 pwm';
saveName = '21Dhh_ASAP2f_ortC1-3_TNT_20pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, ort[C1-3] control (21Dhh, ort-lexA), pwm 20
timeseriesInd = find(i_21Dhh_asap2f_ort...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    (i_21Dhh_asap2f_ort .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort control 20 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort control 20 pwm';
saveName = '21Dhh_ASAP2f_ortC1-3_20pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, TNT control (21Dhh, ort-lexA), pwm 20 
timeseriesInd = find(i_21Dhh_asap2f_TNT...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    (i_21Dhh_asap2f_TNT .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f TNT control 20 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f TNT control 20 pwm';
saveName = '21Dhh_ASAP2f_TNT_20pwm';
inv = 1;
yScale = [-0.1 0.1];


%% L2>>ASAP2f, ort[C1-3]>>TNT (21Dhh, ort-lexA), pwm 200 + 20uM CDM 
% new data from 9/30/17-10/2/17; flyID >= 223 
timeseriesInd = find((flyID>=223).*i_21Dhh_asap2f_ort_TNT_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort TNT + CDM';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort TNT + CDM';
saveName = '21Dhh_ASAP2f_ortC1-3_TNT_200pwm_CDM';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, ort[C1-3] control (21Dhh, ort-lexA), pwm 200 + 20 uM CDM 
% new data from 9/30/17-10/2/17; flyID >= 223 
timeseriesInd = find((flyID>=223).*i_21Dhh_asap2f_ort_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort control + CDM';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort control + CDM';
saveName = '21Dhh_ASAP2f_ortC1-3_200pwm_CDM';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, TNT control (21Dhh, ort-lexA), pwm 200 + 20uM CDM 
% new data from 9/30/17-10/2/17; flyID >= 223 
timeseriesInd = find((flyID>=223).*i_21Dhh_asap2f_TNT_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f TNT control + CDM';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f TNT control + CDM';
saveName = '21Dhh_ASAP2f_TNT_200pwm_CDM';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, ort[C1-3]>>TNT (21Dhh, ort-lexA), pwm 20 + 20uM CDM
% new data from 9/30/17-10/2/17; flyID >= 223 
timeseriesInd = find((flyID>=223).*i_21Dhh_asap2f_ort_TNT_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    ((flyID>=223).*i_21Dhh_asap2f_ort_TNT_CDM .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort TNT 20 pwm + CDM';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort TNT 20 pwm + CDM';
saveName = '21Dhh_ASAP2f_ortC1-3_TNT_20pwm_CDM';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, ort[C1-3] control (21Dhh, ort-lexA), pwm 20 + 20uM CDM 
% new data from 9/30/17-10/2/17; flyID >= 223 
timeseriesInd = find(i_21Dhh_asap2f_ort_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    ((flyID>=223).*i_21Dhh_asap2f_ort_CDM .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f ort control 20 pwm + CDM';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f ort control 20 pwm + CDM';
saveName = '21Dhh_ASAP2f_ortC1-3_20pwm_CDM';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, TNT control (21Dhh, ort-lexA), pwm 20 + 20uM CDM *
% new data from 9/30/17-10/2/17; flyID >= 223 
timeseriesInd = find(i_21Dhh_asap2f_TNT_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    ((flyID>=223).*i_21Dhh_asap2f_TNT_CDM .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f TNT control 20 pwm + CDM';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f TNT control 20 pwm + CDM';
saveName = '21Dhh_ASAP2f_TNT_20pwm_CDM';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, 200pwm (21Dhh)
timeseriesInd = find(i_21Dhh_asap2f...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f 200 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f 200 pwm';
saveName = '21Dhh_ASAP2f_200pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f + 20uM CDM, 200pwm (21Dhh)
timeseriesInd = find(i_21Dhh_asap2f_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f + CDM, 200 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f + CDM 200 pwm';
saveName = '21Dhh_ASAP2f_CDM_200pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f, 20pwm (21Dhh)
timeseriesInd = find(i_21Dhh_asap2f...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    (i_21Dhh_asap2f .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f 20 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f 20 pwm';
saveName = '21Dhh_ASAP2f_20pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f + 20uM CDM, 20pwm (21Dhh)
timeseriesInd = find(i_21Dhh_asap2f_CDM...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    (i_21Dhh_asap2f_CDM .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f + CDM, 20 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f + CDM 20 pwm';
saveName = '21Dhh_ASAP2f_CDM_20pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L1>>ASAP2f, layer M1, 200pwm (GMR37E04)
timeseriesInd = find(i_GMR37E04_asap2f...
    .*~isMoving.*i_M1.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L1 layer M1 ASAP2f, 200 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L1 layer M1 ASAP2f, 200 pwm';
saveName = 'GMR37E04_M1_ASAP2f_200pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L1>>ASAP2f, layer M1, + 20uM CDM, 200pwm (21Dhh)
timeseriesInd = find(i_GMR37E04_asap2f_CDM...
    .*~isMoving.*i_M1.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L1 ASAP2f + CDM, 200 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L1 ASAP2f + CDM 200 pwm';
saveName = 'GMR37E04_M1_ASAP2f_CDM_200pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L1>>ASAP2f, layer M1, 20pwm (21Dhh)
timeseriesInd = find(i_GMR37E04_asap2f...
    .*~isMoving.*i_M1.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    (i_GMR37E04_asap2f .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M1.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L1 layer M1 ASAP2f 20 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L1 layer M1 ASAP2f 20 pwm';
saveName = 'GMR37E04_M1_ASAP2f_20pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L1>>ASAP2f, layer M1, + 20uM CDM, 20pwm (21Dhh)
timeseriesInd = find(i_GMR37E04_asap2f_CDM...
    .*~isMoving.*i_M1.*i_482_18.*(i_ND==1.3).*(i_pwm==20).*...
    i_fullfield_LDflash20ms_Gray500ms + ...
    (i_GMR37E04_asap2f_CDM .* i_300msSearchStimFlash ...
    .*~isMoving.*i_M1.*i_482_18.*(i_ND==1.3).*(i_pwm==200)));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L1 ASAP2f layer M1 + CDM, 20 pwm';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L1 ASAP2f layer M1 + CDM 20 pwm';
saveName = 'GMR37E04_M1_ASAP2f_CDM_20pwm';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f (21Dhh), L2>>TNT (GMR16H03)
timeseriesInd = find(i_21Dhh_asap2f_GMR16H03_TNT...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, TNT';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, TNT';
saveName = '21Dhh_ASAP2f_GMR16H03_TNT';
inv = 1;
yScale = [-0.1 0.1];

%% L2>>ASAP2f (21Dhh), L2 (GMR16H03) control
timeseriesInd = find(i_21Dhh_asap2f_GMR16H03...
    .*~isMoving.*i_M2.*i_482_18.*(i_ND==1.3).*(i_pwm==200).*...
    (i_300msSearchStimFlash + i_fullfield_LDflash20ms_Gray500ms));
refStimCode = '300ms_searchStimFlash';
nrefStimCode = 'fullfield_LDflash20ms_Gray500ms';
pairedEpochs = [1,2]; 
binWidthMult = 1;
interpFrameRate = floor(120); % sample to 120 Hz
refPlotTitle = '300msSearchStimFlash L2 ASAP2f, GMR16H03 control';
nrefPlotTitle = 'LightDarkFlash20ms Gray500ms L2 ASAP2f, GMR16H03 control';
saveName = '21Dhh_ASAP2f_GMR16H03';
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
