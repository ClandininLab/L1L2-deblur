% compute_bootstrappedMetrics.m
%
% Script to load impulse responses from full field flash onto gray saved 
%  data and compute bootstrapped metrics 
%
% Updated: HHY 10/13/17
%

clear all
close all
clc

% number of bootstrap repetitions
nBReps = 10000;
% number of samples to draw per bootstrap
% nSamples = 30;

% end of 2nd phase 0.25 sec from stim start
tEndPhase2 = 0.25;


% savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/CellTypeSpecificQuantificationTest/';
savePath = '/Users/hyang/Documents/New Imaging Analysis Temp/Quantification/';


% get directory
% dirName = uigetdir('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/');
dirName = uigetdir('/Users/hyang/Documents/New Imaging Analysis Temp/');

cd(dirName);

files = dir;

% run processing on every .mat file in directory
for i = 3:length(files)
    clear metricStrctArry t
    
    % load .mat file
    load(files(i).name);
    
    % data from all responding cells
    respROIMat = roiDataMat(iResp,2);
    
    nSamples = length(respROIMat); % number of ROIs
    
    % times corresponding to each frame; only need 1 copy
    t = respROIMat(1).t{1}; 
    ifi = median(diff(t));
    
    % frame corresponding to end of second phase
    endPhase2 = floor(tEndPhase2/ifi);
    
    % get ROI indicies of bootstrapped samples; matrix of nBReps x nSamples
    bootInd = randi(nSamples, nBReps, nSamples);
    
    % compute average traces and metrics for all bootstrapped samples
    for j = 1:nBReps
        metricStrctArry(j).ind = bootInd(j,:); % save indicies used
        
        % compute average traces and save them
        rats1 = zeros(nSamples,length(respROIMat(1).rats{1}));
        rats2 = rats1;
        for n = 1:nSamples
            sampInd = bootInd(j,n);
            rats1(n,:) = respROIMat(sampInd).rats{1}; % dark
            rats2(n,:) = respROIMat(sampInd).rats{2}; % light
        end 
        % avgTrace for dark
        metricStrctArry(j).avgTrace(1,:) = mean(rats1,1); 
        % avgTrace for light
        metricStrctArry(j).avgTrace(2,:) = mean(rats2,1);
        
        % compute metrics
        
        % framePeaks of 1st and 2nd phases for dark
        [metricStrctArry(j).framePeak1(1),...
            metricStrctArry(j).framePeak2(1)]  = computeFramePeaks(...
            metricStrctArry(j).avgTrace(1,:),1);
        % framePeaks for light
        [metricStrctArry(j).framePeak1(2),...
            metricStrctArry(j).framePeak2(2)]  = computeFramePeaks(...
            metricStrctArry(j).avgTrace(2,:),2);
        % convert to column vectors
        metricStrctArry(j).framePeak1 = metricStrctArry(j).framePeak1';
        metricStrctArry(j).framePeak2 = metricStrctArry(j).framePeak2';
        
        % first zero crossing
        % dark
        metricStrctArry(j).frameZero1(1) = computeFrameZero1(...
            metricStrctArry(j).avgTrace(1,:),...
            metricStrctArry(j).framePeak1(1),1);
        % light
        metricStrctArry(j).frameZero1(2) = computeFrameZero1(...
            metricStrctArry(j).avgTrace(2,:),...
            metricStrctArry(j).framePeak1(2),2);
        % convert to column vector
        metricStrctArry(j).frameZero1 = metricStrctArry(j).frameZero1';
        
        % second zero crossing
        % dark
        metricStrctArry(j).frameZero2(1) = computeFrameZero2(...
            metricStrctArry(j).avgTrace(1,:),...
            metricStrctArry(j).framePeak1(1),1);
        % light
        metricStrctArry(j).frameZero2(2) = computeFrameZero2(...
            metricStrctArry(j).avgTrace(2,:),...
            metricStrctArry(j).framePeak1(2),2);
        % convert to column vector
        metricStrctArry(j).frameZero2 = metricStrctArry(j).frameZero2';
        
        % tPeaks - always relative to first zero crossing
        metricStrctArry(j).tPeak1 = t(metricStrctArry(j).framePeak1 - ...
            metricStrctArry(j).frameZero1)';
        metricStrctArry(j).tPeak2 = t(metricStrctArry(j).framePeak2 - ...
            metricStrctArry(j).frameZero1)';
        
        % area under curve of first phase in units %dF/F * sec
        %  between first and second zero crossings
        % dark
        metricStrctArry(j).area1(1) = trapz(...
            metricStrctArry(j).avgTrace(1,...
            metricStrctArry(j).frameZero1(1):...
            metricStrctArry(j).frameZero2(1))) * 100 * ifi;
        % light
        metricStrctArry(j).area1(2) = trapz(...
            metricStrctArry(j).avgTrace(2,...
            metricStrctArry(j).frameZero1(2):...
            metricStrctArry(j).frameZero2(2))) * 100 * ifi;
        % convert to column vector
        metricStrctArry(j).area1 = metricStrctArry(j).area1';
        
        % area under curve of second phase in units %dF/F * sec
        %  between second zero crossing and endPhase2 (0.25sec from stim
        %  start)
        % dark
        metricStrctArry(j).area2(1) = trapz(...
            metricStrctArry(j).avgTrace(1,...
            metricStrctArry(j).frameZero2(1):endPhase2)) * 100 * ifi;        
        % light
        metricStrctArry(j).area2(2) = trapz(...
            metricStrctArry(j).avgTrace(2,...
            metricStrctArry(j).frameZero2(2):endPhase2)) * 100 * ifi;
        % convert to column vector
        metricStrctArry(j).area2 = metricStrctArry(j).area2';
        
        % ratio of 2 area of second phase/area of first phase
        metricStrctArry(j).areaRatio(1) = metricStrctArry(j).area2(1) / ...
            metricStrctArry(j).area1(1);
        metricStrctArry(j).areaRatio(2) = metricStrctArry(j).area2(2) / ...
            metricStrctArry(j).area1(2); 
        % convert to column vector
        metricStrctArry(j).areaRatio = metricStrctArry(j).areaRatio';
    end
    
    % save new .mat file
    save([savePath filesep files(i).name], 'metricStrctArry','t','-v7.3');
end