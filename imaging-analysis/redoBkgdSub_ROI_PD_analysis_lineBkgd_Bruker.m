%% redoBkgdSub_ROI_PD_analysis_lineBkgd_Bruker.m
% 
% Modified script to rerun only background subtraction for data where this
%  has already been done
% No ROI selection. Expects this to already be in pData file
% No processing of stimulus. Expects this to already be in pData file

clc;
clearvars;
close all;

%% Access the time series directory 
disp('Select a time series (TSeries) folder:'); 
seriesDir = uigetdir;
disp(seriesDir);
curDir = pwd;
cd(seriesDir);
m = dir('imDat*');

if (isempty(m)) % exit if the imaging data file does not exist 
    disp('No imaging data (imDat.mat) file found.');
    return;
end 

%% Load previous pData
pDataFile = dir('TSeries*_pData.mat');
load(pDataFile.name);

%% Load imaging data (aligned and unaligned series)
filename = 'imDat.mat';
load(filename ,'unSeries', 'alSeries');

% add these to pDat
pDat.unSeries = unSeries;
pDat.alSeries = alSeries;

%% Calculate background region
[bkMask ] = selectBackground(pDat.avgIm);
pDat.bkMask = bkMask;

%% Get average signal and background-subtracted signal
numLinesPool = 1; % number of lines to pool over 
minPx = 1; % minimum number of pixels to use for estimating background

[avSignal, bksSignal] = getRawSignalsLineBkgdSub(pDat, numLinesPool, minPx);

% add to pDat
pDat.avSignal = avSignal;
pDat.dSignal = bksSignal;

%% Plot signal in ROIs
displayRawSignals(pDat);

%% Save the essentials of pDat and save the stimulus metadata into pDat
% generates a _pData.mat file in the time series directory folder
saveProcessedData_Bruker(pDat, pDat.stimDat);

% go back to prev directory
cd(curDir);
