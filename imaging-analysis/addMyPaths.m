% addMyPaths.m
% 
% script to add relevant 2p data analysis paths
%
% Users: change the following paths to match those on your local computer

%% Bio-Formats package
bioformatsPath = '/Users/hyang/Documents/MATLAB/bfmatlab';
addpath(bioformatsPath);
    
%% Analysis code repository (2p-analysis-code)
analysisPath = '/Users/hyang/Documents/2p-analysis-code';
addpath(analysisPath);

%% More analysis code
hhyAnalysisPath = '/Users/hyang/Documents/2p-analysis-code/HHY_stimulusSpecificAnalysisScripts';
addpath(hhyAnalysisPath);
    
%% Stimulus code repository (2p-stim-code)
stimPath = '/Users/hyang/Documents/2p-stim-code';
addpath(stimPath);
    
%% Folder containing your metadata spreadsheet 
metadataPath = '/Users/hyang/Documents/New Imaging Analysis';
addpath(metadataPath);
    
%% Folder containing processed data (*_pdata.mat files)
pDataPath = '/Users/hyang/Documents/New Imaging Analysis/pData';
addpath(pDataPath);