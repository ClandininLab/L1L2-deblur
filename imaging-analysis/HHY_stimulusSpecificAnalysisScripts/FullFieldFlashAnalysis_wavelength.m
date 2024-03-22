% FullFieldFlashAnalysis_wavelength.m
%
% Full script to analyze full field flash data. Returns data and figures 
%  for a specific genotype.
% This particular variant defines the reference time series by wavelength,
%  not by stimulus
%
% Run after ROI_PD_analysis_singleTimeSeries.m and collect_pdata.m. That
%  is, make sure pData exists for genotypes to be analyzed.
%
% Combines selectTimeSeries.m, stimPlots.m, and
%  analyze_FullFieldFlash_BinaryContrast.m, with a few adaptations. 
%  HHY prefers this organization to splitting it into separate files
%
% Last Updated: 3/13/17
%

close all;
clear all;
clc

% folder for saving ROI data matrix, for stimulus-specific analyses  
figPath = '/Users/hyang/Documents/New Imaging Analysis/Figures/170313';
dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170313';

addMyPaths;

%% Read metadata spreadsheet and find indices 
% struct containing arrays or cell arrays for each column in your 
% experimental metadata spreadsheet 
readDatabase_selectSamples;

%% Select time series of interest 

% REQUIRED VARIABLES (see sections below for examples)
%   timeseriesInd
%   refStimCode (set to '' if you don't have a reference stimulus)
%   refPlotTitle
%   nrefPlotTitle (e.g. title for FFFOntoGray plots)
%   yScale 
%   inv  - set to 1 for ASAP voltage imaging
%   pairedEpochs - for FFFOntoGray only to indicate light/dark epoch pairs
%   saveName - if you want the code to automatically save your data matrices


%% L2>>jRGECO1b, Ace2N-mNeon
timeseriesInd = find(i_21Dhh_ace2nmneon_jrgeco1b.*i_M2.*i_300msFFF...
    .*(~metaDat.isMoving).*(metaDat.estFrameRate==40));
refWavelength = 1050;
nRefWavelength = [920 950 980 1010];
refPlotTitle = 'L2 jRGECO1b';
plotTitle = 'L2 Ace2N-mNeon';
interpFrameRate = floor(120); % sample to 120 Hz
binWidthMult = 3; % sliding average
saveName = 'L2_Ace2NmNeon_jRGECO1b_300msFFF';
inv = 0; 
refYScale = [-0.3 0.3];
nRefYScale = [-0.015 0.015];

%% Load ROI data, doing ROI-matching across time series if necessary

% If no reference stim, create an R x S matrix where R is the number ...
% of ROIs and S is the number of stimuli per ROI
if isempty(refWavelength)
    roiDataMat = createROIMatrix(timeseriesInd, metaDat, pDataPath);
    
else % reference stim exists 
    % match ROIs across different time series
    roiMetaMat = matchROIsAcrossWavelength(timeseriesInd, refWavelength, ...
        metaDat);
    
    % % view struct array field values across ROIs
    % {roiMetaMat(:,1).stimclass}';

    % for each ROI, load the pData and save the experimental metadata into the
    % struct. 
    roiDataMat = loadROIData(roiMetaMat, metaDat, pDataPath);
end 

% clear roiMetaMat metaDat

%% Analyze FullFieldFlash_BinaryContrast ref stimulus responses 
% refColumn = 1;
% % if roiDataMat has only one column, then the 'reference' column is 1
% if size(roiDataMat,2) ~= 1
%     roiDataMatTemp = roiDataMat;
%     for r = 1:size(roiDataMat,1)
%         for s = 1:size(roiDataMat,2)
%             currentWavelength = roiDataMatTemp(r,s).wavelength;
%             if (currentWavelength==refWavelength)
%                 % ref column is always the first column
%                 roiDataMat(r,1) = roiDataMatTemp(r,s);
%             else
%                 % rest of columns are in order defined by nRefWavelength
%                 currentInd = find(nRefWavelength == currentWavelength);
%                 roiDataMat(r,currentInd+1) = roiDataMatTemp(r,s);
%             end 
%         end 
%     end
% end

[refColumn, roiDataMat] = sortROIDataMatByWavelength(roiDataMat,...
    refWavelength,nRefWavelength);

% screen responses of individual ROIs, label as responding,
% non-responding, or inverted
[roiDataMat, iRefResp, iRefInv, BIN_WIDTH] = ...
    aggregate_fff_matchedROIs(roiDataMat, refColumn, interpFrameRate, ...
    binWidthMult);

% filter out all ROIs in the data matrix based on reference ROIs 
respROIMat = roiDataMat(iRefResp,:);


%% Plots
plot_FullFieldFlash(respROIMat(:,1),inv,refYScale,1/interpFrameRate,...
    roiDataMat(1,1).stimDat.Duration{1},...
    [refPlotTitle ' ' num2str(refWavelength) 'nm'],0);
for i = 2:size(respROIMat,2)
    flashDuration = roiDataMat(1,i).stimDat.Duration{1};
    plot_FullFieldFlash(respROIMat(:,i), inv, nRefYScale, ...
        1/interpFrameRate, flashDuration, ...
        [plotTitle ' ' num2str(nRefWavelength(i-1)) 'nm'], 0);
end

%% Save stimulus-specific params
save([dataPath filesep saveName], 'respROIMat', 'roiDataMat', ...
    'iRefResp', 'iRefInv', 'inv', 'refPlotTitle', ...
    'refYScale','nRefYScale','plotTitle','nRefWavelength');
% saveas(figID_idv,[figPath filesep saveName '_idv'],'fig');
% saveas(figID_avg,[figPath filesep saveName '_avg'],'fig');
