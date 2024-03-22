% convert3ContrastFileFormat.m
%
% Script to take 3 contrast responses from previous imaging setup and
% convert it to a file format that can be run through the new code base, to
% compute bootstrapped metrics and to plot
%
% Updated: 11/22/17 HHY
%

clearvars
close all
clc

lightInd = 1:63;
darkInd = 64:126;

dataFolder = '/Users/hyang/Documents/New Imaging Analysis Temp/';
cd(dataFolder)

% get old data file to convert
[fileName, pathName, ~] = uigetfile('.mat');
load([pathName fileName]);
% remove .mat from file name
newFileName = fileName(1:(end-4));

% loop through 3 different contrasts
for i = 1:length(ROIresponsesShift.contrastStimCodes)
    % define which contrast
    switch i
        case 1
            contrastCode = '05c';
        case 2
            contrastCode = '025c';
        case 3
            contrastCode = '0125c';
    end
    
    iResp = [];
    roiDataMat = [];
    
    % loop through individual cells, reassign to new data format
    for j = 1:length(ROIresponsesShift.responses)
        % responding?
        if (ROIresponsesShift.responses(j).responding == 1)
            iResp = [iResp j];
        end
        
        % roiDataMat for refStim
        roiDataMat(j,1).rats = ...
            ROIresponsesShift.responses(j).refStim.rates';
        roiDataMat(j,1).stdErr = ...
            ROIresponsesShift.responses(j).refStim.stdErrs';
        roiDataMat(j,1).stimcode = ...
            ROIresponsesShift.refStimCode;
        roiDataMat(j,1).flyID = ...
            ROIresponsesShift.responses(j).flyID;
        
        % roiDataMat for short flash stim
        % dark = 1, light = 2
        roiDataMat(j,2).rats{1} = ...
            ROIresponsesShift.responses(j).contrastStim(i).rates(darkInd);
        roiDataMat(j,2).rats{2} = ...
            ROIresponsesShift.responses(j).contrastStim(i).rates(lightInd);
        roiDataMat(j,2).stdErrs{1} = ...
            ROIresponsesShift.responses(j).contrastStim(i).stdErrs(darkInd);
        roiDataMat(j,2).stdErrs{2} = ...
            ROIresponsesShift.responses(j).contrastStim(i).stdErrs(lightInd);
        roiDataMat(j,2).flyID = ...
            ROIresponsesShift.responses(j).flyID;
        roiDataMat(j,2).seriesID = ...
            ROIresponsesShift.responses(j).contrastStim(i).name;
        roiDataMat(j,2).stimcode = ...
            ROIresponsesShift.contrastStimCodes{i};
        
        roiDataMat(j,2).stimDat.FlashDuration = ...
            ROIresponsesShift.flashTime;
        roiDataMat(j,2).stimDat.GrayDuration = ...
            ROIresponsesShift.interTime;
        
        tempT = 0:(length(roiDataMat(j,2).rats{1})-1);
        roiDataMat(j,2).t{1} = tempT * intIFI;
        tempT = 0:(length(roiDataMat(j,2).rats{2})-1);
        roiDataMat(j,2).t{2} = tempT * intIFI;
    end

    % save data file just for that contrast
    
    save([pathName newFileName '_' contrastCode], 'roiDataMat', 'iResp', '-v7.3');
    
end


