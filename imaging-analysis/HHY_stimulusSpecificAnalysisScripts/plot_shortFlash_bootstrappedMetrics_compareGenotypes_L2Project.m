% plot_shortFlash_bootstrappedMetrics_compareGenotypes_L2Project.m
%
% Script to generate plots of bootstrapped metrics (from
%  compute_bootstrappedMetrics.m) for light and dark impulse responses.
%  Quantification for the L1/L2 project. Plots mean and specified
%  confidence interval.
%
% Updated: HHY 10/13/17
%
%

clearvars
close all
clc

% cd ('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/CellTypeSpecificQuantification');
% cd ('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/ContrastDataQuantification');
cd ('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/Quantification2 - this one in paper');

ci = 95;
% ci = 99;
% ci = 99.9;

ciLowVal = ((100 - ci)/2)/100;
ciHighVal = 1 - ciLowVal;

numGeno = 3;

fieldList = {'area1', 'area2', 'areaRatio', 'tPeak1', 'tPeak2'};
% fieldList = {'tPeak2'};


yScale = [0 0.15; -0.05 0.1; -1 2; 0 0.05; 0 0.55];
% yScale = [0 0.2];

%% Load Data
% preallocate
dataMatrix = zeros(10000,numGeno);

for k=1:numGeno
    [filename, pathname] = uigetfile('*.mat');

    % load data
    load([pathname filename], 'metricStrctArry');
    
    msa{k} = metricStrctArry;
    
    names{k} = filename;
    
    clear metricStrctArry
end
%% Plot all Figures

% light and dark
for i=1:2
    % plot all metrics
    for j=1:length(fieldList)
        % get all genotypes
        for k=1:numGeno
            % temp save metric
            dataMatrix(:,k) = returnMetricSAF(msa{k},fieldList{j},i);
        end
        % invert some metrics
        if ((i==1)&&(strcmpi(fieldList{j},'area1')))||...
                ((i==2)&&(strcmpi(fieldList{j},'area2')))||...
                (strcmpi(fieldList{j},'areaRatio'))
            dataMatrix = dataMatrix * -1;
        end
        
        % compute mean
        allMeans = mean(dataMatrix,1);
        % compute confidence interval
        ciLowAll = allMeans - quantile(dataMatrix,ciLowVal,1);
        ciHighAll = quantile(dataMatrix,ciHighVal,1) - allMeans;
        
        % make figure
        figure; hold on;
        for k=1:numGeno
            h(k) = errorbar(k,allMeans(k),ciLowAll(k),ciHighAll(k),'o');
        end
        if (i==1)
            contrast = 'dark';
        elseif (i==2)
            contrast = 'light';
        end
        title(['ci= ' num2str(ci) ' ' fieldList{j} ' ' contrast]);
        legend(h,names);
        
        ylim(yScale(j,:));
    end
end

%% Plot comparison between light and dark of 1 genotype
% Load Data
[filename, pathname] = uigetfile('*.mat');
load([pathname filename], 'metricStrctArry');

% plot all metrics
for j=1:length(fieldList)
    for i=1:2 % light and dark
        % invert some metrics
        if ((i==1)&&(strcmpi(fieldList{j},'area1')))||...
                ((i==2)&&(strcmpi(fieldList{j},'area2')))||...
                (strcmpi(fieldList{j},'areaRatio'))
            dataMatrix(:,i) = -1*...
                returnMetricSAF(metricStrctArry,fieldList{j},i);
        else
            dataMatrix(:,i) = returnMetricSAF(metricStrctArry,...
                fieldList{j},i);
        end 
    end

    % compute mean
    allMeans = mean(dataMatrix,1);
    % compute confidence interval
    ciLowAll = allMeans - quantile(dataMatrix,ciLowVal,1);
    ciHighAll = quantile(dataMatrix,ciHighVal,1) - allMeans;

    % make figure
    figure; hold on;
    for k=1:2
        h(k) = errorbar(k,allMeans(k),ciLowAll(k),ciHighAll(k),'o');
    end
    title(['ci= ' num2str(ci) ' ' fieldList{j}]);
    legend(h,'Dark', 'Light');

    ylim(yScale(j,:));
end