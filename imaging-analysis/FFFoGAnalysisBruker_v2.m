% FFFoGAnalysisBruker_v2.m
%
% Script for analyzing impulse responses, with option to pool
%
% USER NEEDS TO EDIT EVERYTHING BETWEEN:
%   "(3) Select time series of interest"
%   and
%   "(4) Load ROI data"
%
% also make sure correct version of readDatabase_selectSamples_NatStim and 
%   generateSaveNames are being used
%
% v2: 220607 MMP added option to skip resampling/interpolation, using 
%   filterROIs_simple_stimLockedPlots_FFFoG_pool (did not write analogous
%   function for situation w/ no reference stimulus, need to add if desired) 
% also added if statement to resolve error with sortROIDataMatByStim

%% (1) Set up the workspace

close all;
clearvars;

% pData directory
addMyPaths;
dataPath = savePath;
pdataPath = pDataPath;

%% Read metadata spreadsheet and find indices 
% struct containing arrays or cell arrays for each column in your 
% experimental metadata spreadsheet 

%readDatabase_selectSamples_NatStim;
readDatabase_selectSamples_NatStim_MMP;
who

%% (3) Select time series of interest - EDIT THIS SECTION!

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

%% !!!EDIT THESE TO SELECT TIME SERIES BY ATTRIBUTES!!!

i_flyID = (flyID > 14);
%i_flyID = (flyID == 1);

timeseriesInd = find(...
    (i_21Dhh_ASAP2f_jRGECO1b_new + i_21Dhh_ASAP2f_jRGECO1b_newer) .* ...
    i_flyID .* ...
(i_300msSearchStimFlash_RGB + i_fullfield_LDflash20ms_Gray500ms_RGB))

%    i_M1 .*...
%    (~metaDat.isMoving) .*...

interpFrameRate = -1; % set to -1 if you don't want to resample
yScale = [-0.1 0.1];
inv = 1;

% how many trials of nref stim presented for FOV
trialsPerROI = 1; 

%% may need to edit, may not
refStimCode = '300ms_searchStimFlash_RGB';
%nrefStimCode = {'300msFFF_binarycontrast_RGB'};
refColumn = 1;
nrefColumn = 2;
binWidthMult = 1;
pairedEpochs = [1,2]; 

%% auto-generate save names (220518 MMP)

generateSaveNames;
% refPlotTitle = '300ms_searchStimFlash_RGB';
% nrefPlotTitle = '300msFFF_binarycontrast_RGB';
% saveName = '300msFFF_binarycontrast_RGB';

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
    % throwing an error, so skipping via if statement b/c usually not necessary...
    if length(nrefStimCode)>1
        [refColumn, roiDataMat] = sortROIDataMatByStim(roiDataMat,...
            refStimCode,nrefStimCode);
    end
end 

%% Save data matrix (optional)
save([savePath '/' saveName], 'roiDataMat');

%% (1) Get mean stimulus-locked response for each ROI, select responding ROIs
% pooled across multiple trials

if (interpFrameRate == -1) % -1 means do not interpolate/resample (MMP 220607)
% NEED TO ADD SIMPLIFIED VERSION OF NO-REF IF WE WANT TO USE IT
%     if isempty(refStimCode) % if ref doesn't exist
%         [roiDataMatMeans, iResp, iInv, framesPerLDCycle] = ...
%             filterROIs_stimLockedPlots_FFFoG_noRef_pool(roiDataMat,pairedEpochs,...
%             interpFrameRate, inv, yScale, binWidthMult);
%     else % if reference exists
        [roiDataMatMeans, iResp, iInv, framesPerLDCycle] = ...
            filterROIs_simple_stimLockedPlots_FFFoG_pool(roiDataMat, refColumn, ...
            pairedEpochs, interpFrameRate, inv, yScale, binWidthMult);
%     end

else % original interpolation/resampling pipeline
    if isempty(refStimCode) % if ref doesn't exist
        [roiDataMatMeans, iResp, iInv, framesPerLDCycle] = ...
            filterROIs_stimLockedPlots_FFFoG_noRef_pool(roiDataMat,pairedEpochs,...
            interpFrameRate, inv, yScale, binWidthMult);
    else % if reference exists
        [roiDataMatMeans, iResp, iInv, framesPerLDCycle] = ...
            filterROIs_stimLockedPlots_FFFoG_pool(roiDataMat, refColumn, ...
            pairedEpochs, interpFrameRate, inv, yScale, binWidthMult);
    end 
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
% all cells
plot_FFFOntoGray(roiDataMatMeans(:,nrefColumn), inv, yScale, nrefPlotTitle, ...
    pairedEpochs);
% only responding cells
plot_FFFOntoGray(respROIMat(:,nrefColumn), inv, yScale, nrefPlotTitle, ...
    pairedEpochs);
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
