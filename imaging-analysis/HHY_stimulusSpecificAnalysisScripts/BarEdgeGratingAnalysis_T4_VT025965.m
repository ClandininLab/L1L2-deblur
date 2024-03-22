% BarEdgeGratingAnalysis_v1.m
%
% Script for analyzing responses to moving bar, moving edge, static and
%  moving grating (e.g. for T4/T5 imaging)
% combines selectTimeSeries with stimulus-specific analyses
%
% 6/6/17 - quick analysis of T4 in VT025965

%% (1) Set up the workspace

close all;
clear all;

addMyPaths;

savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170606';
figPath = '/Users/hyang/Documents/New Imaging Analysis/Figures/170606';
filename = '/Users/hyang/Documents/New Imaging Analysis/2p_imaging_log.txt';

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

% %% T4/T5 para FlpStop - experimental
% timeseriesInd = find(i_R59E08AD_R42F06DBD_gcamp6f_paraFlpStop_Flp...
%     .*~isMoving.*i_male.*...
%     (i_movingbar_dark_c1_5deg_4dir_20dps + i_movingedge_4dir_20dps + ...
%     i_staticmovinggrating_20deg_20dps_c1_12dir));
% refStimCode = 'movingBar_darkC1_5deg_4dir_20ds';
% nrefStimCode = {'movingEdge_4dir_20ds', ...
%     'staticMovingGrating_20deg_20ds_c1_12dir'};
% binWidthMult = 1;
% interpFrameRate = floor(20); % sample to 20 Hz
% plotTitle = 'T4/T5 para FlpStop experimental';
% saveName = 'T4T5_gcamp6f_paraFlpStop_exp';
% inv = 0;
% yScale = [-0.5 2.5];
% 
% %% T4/T5 para FlpStop - no Flp control
% timeseriesInd = find(i_R59E08AD_R42F06DBD_gcamp6f_paraFlpStop...
%     .*~isMoving.*i_male.*...
%     (i_movingbar_dark_c1_5deg_4dir_20dps + i_movingedge_4dir_20dps + ...
%     i_staticmovinggrating_20deg_20dps_c1_12dir));
% refStimCode = 'movingBar_darkC1_5deg_4dir_20ds';
% nrefStimCode = {'movingEdge_4dir_20ds', ...
%     'staticMovingGrating_20deg_20ds_c1_12dir'};
% binWidthMult = 1;
% interpFrameRate = floor(20); % sample to 20 Hz
% plotTitle = 'T4/T5 para FlpStop no Flp control';
% saveName = 'T4T5_gcamp6f_paraFlpStop_noFlpCtrl';
% inv = 0;
% yScale = [-0.5 2.5];
% 
% %% T4/T5 para FlpStop - heterozygous female control
% timeseriesInd = find(i_R59E08AD_R42F06DBD_gcamp6f_paraFlpStop_Flp_het...
%     .*~isMoving.*i_female.*...
%     (i_movingbar_dark_c1_5deg_4dir_20dps + i_movingedge_4dir_20dps + ...
%     i_staticmovinggrating_20deg_20dps_c1_12dir));
% refStimCode = 'movingBar_darkC1_5deg_4dir_20ds';
% nrefStimCode = {'movingEdge_4dir_20ds', ...
%     'staticMovingGrating_20deg_20ds_c1_12dir'};
% binWidthMult = 1;
% interpFrameRate = floor(20); % sample to 20 Hz
% plotTitle = 'T4/T5 para FlpStop het control';
% saveName = 'T4T5_gcamp6f_paraFlpStop_hetCtrl';
% inv = 0;
% yScale = [-0.5 2.5];

%% T4 in M10 in VT025965 pattern
timeseriesInd = find(i_VT025965_gcamp6f.*~isMoving.*i_M10.*...
    i_movingbar_light_c1_5deg_4dir_20dps);
refStimCode = '';
nrefStimCode = '';
binWidthMult = 1;
interpFrameRate = floor(40); % sample to 40 Hz
plotTitle = 'T4 dendrites in VT025965';
saveName = 'T4-M10_VT025965_gcamp6f_movingLightBar';
inv = 0;
yScale = [-0.5 2.5];
epochOrder = [1 3 2 4];

%% (4) Load ROI data, doing ROI-matching across time series if necessary

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

% clear roiMetaMat metaDat

% %% Save data matrix (optional)
% save([savePath filesep saveName], 'roiDataMat');

%% Screen responses of individual ROIs, label as responding or not
% Moving Bar/Moving Edge
[roiDataMat, iResp, iInv, iLayer] = ...
    filterROIs_movingBar_noRef(roiDataMat,...
    interpFrameRate, binWidthMult);

respDataMat = roiDataMat(iResp);
    
%% circshift to align all responses - moving bar   
respDataMat = alignResp_movingBar(respDataMat,0.5);

%% circshift to align all responses- moving edges
respDataMat = alignResp_movingBar(respDataMat,0.65);

%% Screen responses of individual ROIs, label as responding or not
% Static and Moving gratings
[roiDataMat, iResp, iInv, iLayer] = ...
    filterROIs_staticMovingGrating_noRef(roiDataMat,...
    interpFrameRate, binWidthMult);

respDataMat = roiDataMat(iResp);

% compute indicies for polar plot (peak dF/F during static, moving)

respDataMat = computePolarInd_staticMovingGrating(respDataMat,...
    interpFrameRate);
    
%% Plot mean responses - moving bar/moving edge

plot_movingBar(respDataMat, iLayer, yScale, plotTitle, epochOrder, ...
    interpFrameRate, [6.5 6.5 6.5 6.5]);

%% Plots - static and moving gratings - mean responses, polar plots

plot_gratings_polar(respDataMat, iLayer, yScale, plotTitle);

%% Save data
save([savePath filesep saveName], 'roiDataMat', 'iResp', 'interpFrameRate',...
    'iLayer','respDataMat','binWidthMult','inv','yScale', '-v7.3');