% searchProcessed_saveMean.m
%
% Script to load processed 300 ms search stim flash data and to save just 
%  the mean response.
% Only for data format for L1/L2 paper.
% To send data to Shaul regarding request by email on 3/10/18.
% Updated to save individual responses as well as the mean.
%
% Updated: 3/14/18

% where to save new .mat files
saveDir = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/SearchStimResponsesForShaul/';

fileDir = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/CurrentAnalyzedData';
cd(fileDir)

% loop through every file in folder
fileList = dir('*.mat');

for i = 1:length(fileList)
    
    % get short flash data file
    fileName = fileList(i).name;
    
    % load data file - only roiDataMat and iResp variables needed
    load(fileName, 'roiDataMat', 'iResp');
    % short flash stimulus, responding cells only
    respROIMat = roiDataMat(iResp, 1);
    
    % preallocate
    rats = zeros(length(respROIMat),...
        length(respROIMat(1).rats));

    % get all individual ROI responses
    for n = 1:length(respROIMat)
        rats(n,:) = respROIMat(n).rats';
    end 
    % mean response
    meanResp = mean(rats,1);
    % individual responses
    indivResp = rats;
    
    % time
    t = 0:(length(respROIMat(1).rats)-1);
    t = t * 1/120; % 120 Hz; why isn't there a t for this stim
    
    % save data file
    save([saveDir fileName], 'meanResp','indivResp','t');
end
    