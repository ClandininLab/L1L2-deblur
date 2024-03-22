%% ROI_PD_analysis_singleTimeSeries.m
% 
% For a single time series, allows user to select ROIs and background
%  threshold. Calculates raw signal (mean pixel intensity) in each ROI and
%  performs background subtraction. Extracts timing information about the
%  imaging and the stimulus from the photodiode signal.
% Saves essential imaging and stimulus data into a new struct, pDat, 
% for storing pre-processed data. pDat is saved in a pData.mat file for the
% respective time series. 
%
% This is a script. Running the whole thing through will do all of the
%  above. Run this after running preprocess.
% 
% Previous versions: 
%   ROI_analysis_singleTimeSeries.m
%   today_single_1ch_v3.m
% 
% Updates
%   3/7/17 - some updates to functions called by this script, only one that
%       changes what this script sees is the removal of the output of
%       displayRawSignals
%   3/8/17 - more comments
%   3/28/17 - slight edit to processPD
% Last updated: 3/28/17 HHY
%

clc;
clear all;
close all;

%% Access the time series directory 
display('Select a time series (LDM) folder:'); 
seriesDir = uigetdir;
disp(seriesDir);
curDir = pwd;
cd(seriesDir);
m = dir('imDat*');

if (isempty(m)) % exit if the imaging data file does not exist 
    display('No imaging data (imDat.mat) file found.');
    return;
end 

%% Load imaging metadata
filename = 'imDat.mat';
pDat = load(filename);
% fields = {'fileloc', 'name', 'avgIm', 'nFrames', 'xPixels', 'yPixels'};
% pDat = load(filename, fields{:});

%% Play a movie of the unaligned or the aligned time series to evaluate  
% degree movement in this series
% display('Watch a movie of the aligned time series. Mark this series as moving if needed'); 
%     fig = figure; 
%      for f = 1:totFrames
%           imagesc(unalignedSeries(:,:,f)); colormap('gray'); 
%           pause(0.05);
%      end 
%      close(fig);

%% Draw ROIs using the average image of the time series
[roiMasks, numMasks] = drawROIs(pDat.avgIm);
pDat.roiMasks = roiMasks;
pDat.nMasks = numMasks;  

%% Calculate background region
[bkMask ] = selectBackground(pDat.avgIm);
pDat.bkMask = bkMask;

%% Get average signal and background-subtracted signal
[avSignal, bksSignal] = getRawSignals(pDat);
pDat.avSignal = avSignal;
pDat.dSignal = bksSignal;
    
%%  Load stimulus metadata 
stimDat = load('stim.mat');
% extract imaging timing information 
[ifi, imFrameStartTimes] = getImagingIFI(stimDat.obj, 0); 
pDat.imFrameStartTimes = imFrameStartTimes(1:length(pDat.dSignal));
pDat.imIFI = ifi;

%% Process photodiode signals 
[stimDat, stimvalIF, stimEpochTimesIF, rcStim, rcStimInd, ...
    stimEpochStartTimes, lightStartTimes, darkStartTimes, pdThresh] ...
    = processPD(stimDat, pDat.imFrameStartTimes);

% save output into stimulus info struct
pStimDat.stimvalIF = stimvalIF;
pStimDat.stimEpochTimesIF = stimEpochTimesIF;
pStimDat.rcStim = rcStim;
pStimDat.rcStimInd = rcStimInd; 
pStimDat.stimEpochStartTimes = stimEpochStartTimes;
pStimDat.lightStartTimes = lightStartTimes;
pStimDat.darkStartTimes = darkStartTimes;
pStimDat.pdThresh = pdThresh;

% save this into main pDat struct
pDat.pStimDat = pStimDat;

%% Plot signal in ROIs
displayRawSignals(pDat);

%% Save the essentials of pDat and save the stimulus metadata into pDat
% generates a _pData.mat file in the time series directory folder
saveProcessedData(pDat, stimDat);




