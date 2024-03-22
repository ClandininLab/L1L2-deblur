% compute_bootstrappedMetrics_whiteNoise.m
%
% Script to compute bootstrapped metrics (as in
% compute_bootstrappedMetrics.m) only on white noise or flash responses
% averaged to approximate
%
% Updated: HHY 10/13/17

clear all
close all
clc

% number of bootstrap repetitions
nBReps = 10000;

% end of 2nd phase 0.25 sec from stim start
tEndPhase2 = 0.25;

savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/Quantification2/';

%% white noise
cd('/Users/hyang/Documents/New Imaging Analysis/Marjorie/AnalyzedData/whiteNoise');
filename = 'fff300ms_fffrand9_16ms_L2ASAP2f_refFiltered_filtered_analyzed_170915.mat';
load(filename);
k = 44;
linFiltersPlot = -1*linFilters(1:k,useROIs);
t = (ifi*-0.52):ifi:((k-1.52)*ifi);
responses = linFiltersPlot';
m1 = mean(linFiltersPlot,2);

%% short flashes
cd('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170630');
filename = '21Dhh_ASAP2f_200pwm.mat';
load(filename);

respROIMat = roiDataMat(iResp,2);

meanIR = zeros(length(respROIMat), length(respROIMat(1).rats{1}));
for n = 1:length(respROIMat)
    meanIR(n, :) = (respROIMat(n).rats{1} - respROIMat(n).rats{2})/2;
end

m2 = mean(meanIR,1);
responses = meanIR / max(abs(m2)) * max(abs(m1));
t = respROIMat(1).t{1};

%% bootstrap
ifi = median(diff(t));

nSamples = size(responses,1); % number of ROIs

% frame corresponding to end of second phase
endPhase2 = floor(tEndPhase2/ifi);

% get ROI indicies of bootstrapped samples; matrix of nBReps x nSamples
bootInd = randi(nSamples, nBReps, nSamples);

% compute average traces and metrics for all bootstrapped samples
for j = 1:nBReps
    metricStrctArry(j).ind = bootInd(j,:); % save indicies used

    % compute average traces and save them
    rats = zeros(nSamples,size(responses,2));
    for n = 1:nSamples
        sampInd = bootInd(j,n);
        rats(n,:) = responses(sampInd,:);
    end 
    % avgTrace
    metricStrctArry(j).avgTrace = mean(rats,1); 

    % compute metrics

    % framePeaks of 1st and 2nd phases
    [metricStrctArry(j).framePeak1,...
        metricStrctArry(j).framePeak2]  = computeFramePeaks(...
        metricStrctArry(j).avgTrace,1);

    % first zero crossing
    metricStrctArry(j).frameZero1 = computeFrameZero1(...
        metricStrctArry(j).avgTrace,...
        metricStrctArry(j).framePeak1,1);

    % second zero crossing
    metricStrctArry(j).frameZero2 = computeFrameZero2(...
        metricStrctArry(j).avgTrace,...
        metricStrctArry(j).framePeak1,1);

    % tPeaks - always relative to first zero crossing
    metricStrctArry(j).tPeak1 = t(metricStrctArry(j).framePeak1 - ...
        metricStrctArry(j).frameZero1);
    metricStrctArry(j).tPeak2 = t(metricStrctArry(j).framePeak2 - ...
        metricStrctArry(j).frameZero1);

    % area under curve of first phase in units %dF/F * sec
    %  between first and second zero crossings
    metricStrctArry(j).area1 = trapz(...
        metricStrctArry(j).avgTrace(...
        metricStrctArry(j).frameZero1:...
        metricStrctArry(j).frameZero2)) * 100 * ifi;

    % area under curve of second phase in units %dF/F * sec
    %  between second zero crossing and endPhase2 (0.25sec from stim
    %  start)
    metricStrctArry(j).area2 = trapz(...
        metricStrctArry(j).avgTrace(...
        metricStrctArry(j).frameZero2:endPhase2)) * 100 * ifi;

    % ratio of 2 area of second phase/area of first phase
    metricStrctArry(j).areaRatio = metricStrctArry(j).area2 / ...
        metricStrctArry(j).area1; 
end

% save new .mat file
save([savePath filesep filename], 'metricStrctArry','t','-v7.3');