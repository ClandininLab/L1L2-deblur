% computeSignificance_shortFlash_bootstrappedMetrics_L2Project.m
%
% Script to print to screen significance cutoffs of bootstrapped metrics 
%  (from compute_bootstrappedMetrics.m) for light and dark impulse 
%  responses. Pairwise comparisons for selected genotypes. 
%  Quantification for the L1/L2 project. Plots mean and specified
%  confidence interval.
% Note: currently only works for 2 genotypes selected. 
%
% Updated: HHY 10/13/17
%

clearvars
close all
clc

% cd ('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/CellTypeSpecificQuantificationTest');
cd ('/Users/hyang/Documents/New Imaging Analysis Temp/Quantification');

numGeno = 2;
% ci = [95 99 99.9];
ci = [97.5 99.5 99.95]; % for multiple comparisons

fieldList = {'area1', 'area2', 'areaRatio', 'tPeak1', 'tPeak2'};

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

%% Print significance cutoffs
% display which pairs compared
fprintf('%s vs. %s\n', names{1}, names{2});
for i=1:2
    % display light or dark
    switch i
        case 1
            fprintf('Dark\n');  
        case 2
            fprintf('Light\n');
    end
    
    for j=1:length(fieldList)
        % display field
        fprintf('%s\n',fieldList{j});
        
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
        
        % compute means
        allMeans = mean(dataMatrix,1);
        
        for l = 1:length(ci)
            % compute confidence intervals
            ciLowVal = ((100 - ci(l))/2)/100;
            ciHighVal = 1 - ciLowVal;
            ciLowAll = quantile(dataMatrix,ciLowVal,1);
            ciHighAll = quantile(dataMatrix,ciHighVal,1);
            
            % so far, only works for pairwise comparisons
            if (ciLowAll(1) > ciHighAll(2))||...
                    (ciHighAll(1) < ciLowAll(2))
                fprintf('CI %2.1f sig\n',ci(l));
            else
                fprintf('CI %2.1f no sig\n',ci(l));
            end
        end 
    end
end

%% For light vs dark of 1 genotype
% Load Data
[filename, pathname] = uigetfile('*.mat');
load([pathname filename], 'metricStrctArry');

% Print significance cutoffs

% display file name
fprintf('%s Light vs Dark\n',filename);

% preallocate
dataMatrix = zeros(10000,2);
    
for j=1:length(fieldList)
    % display field
    fprintf('%s\n',fieldList{j});

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



    % compute means
    allMeans = mean(dataMatrix,1);

    for l = 1:length(ci)
        % compute confidence intervals
        ciLowVal = ((100 - ci(l))/2)/100;
        ciHighVal = 1 - ciLowVal;
        ciLowAll = quantile(dataMatrix,ciLowVal,1);
        ciHighAll = quantile(dataMatrix,ciHighVal,1);

        % so far, only works for pairwise comparisons
        if (ciLowAll(1) > ciHighAll(2))||...
                (ciHighAll(1) < ciLowAll(2))
            fprintf('CI %2.1f sig\n',ci(l));
        else
            fprintf('CI %2.1f no sig\n',ci(l));
        end
    end 
end