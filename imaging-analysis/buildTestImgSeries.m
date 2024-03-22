% to build test image series for getRawSignalsLineBkgdSub.m
% when writing that code
%
% CREATED: 8/4/22 - HHY

numLines = 50;
numPxPerLine = 200;
nFrames = 3;

% intensity values
roiIntVal(1) = 500;
roiIntVal(2) = 300;
roiIntVal(3) = 10;
noiseVal(1) = 20;
noiseVal(2) = 100;
noiseVal(3) = 0;

% noise stripe pattern, by index into noiseVal
noisePattern = randi(3,numLines,nFrames);
% roi pattern
% roiPattern = randi(3,nFrames,1);
roiPattern = [1 2 3];


quarterImg = zeros(numLines, numPxPerLine / 4);

% loop through each line of quarterImg, with each line, flip 1 more px to
%  1
for i = 1:numLines
    quarterImg(i,1:i) = 1;
end

% use quarter image to build up full image: black 1/8, quarterImg, black
% 1/4, quarterImg, black 1/8
fullImg = [zeros(numLines, numPxPerLine / 8) quarterImg ...
    zeros(numLines, numPxPerLine / 4) quarterImg ...
    zeros(numLines, numPxPerLine / 8)];

fullImgSeries = zeros(numLines, numPxPerLine, nFrames);
expectSignal = zeros(nFrames,1);

% convert fullImg to image series
for i = 1:nFrames
    fullImgSeries(:,:,i) = fullImg * roiIntVal(roiPattern(i));
    expectSignal(i) = roiIntVal(roiPattern(i));
end

% ROI masks - perfectly covers "cells"
roiLog = fullImg > 0;

maskLeft = false(numLines, numPxPerLine);
maskLeft(:,1:(numPxPerLine/2)) = true;
maskRight = false(numLines, numPxPerLine);
maskRight(:,(numPxPerLine/2 + 1):end) = true;

in.roiMasks{1} = roiLog & maskLeft;
in.roiMasks{2} = roiLog & maskRight;

% add stripey noise to image series
fullImgSeriesStripes = zeros(size(fullImgSeries));

% loop through all frames
for i = 1:nFrames
    % loop through all lines
    for j = 1:numLines
        fullImgSeriesStripes(j,:,i) = fullImgSeries(j,:,i)+ ...
            noiseVal(noisePattern(j,i));
    end
end

% background mask
bkgdMask = false(numLines,numPxPerLine);

% BkgdMask1: manually defined, 2 columns, 1 on each side, 10 px wide
%  first column covers lines 1 to 30, second column lines 25 to 5-
% bkgdMask(1:30,1:10) = true;
% bkgdMask(25:50,190:200) = true;

% BkgdMask2: manually defined 2 columns, 1 on each side
%  first column, lines 1 to 40, columns 1 to 10, 10 pixels wide,
%  second column, lines 30 to 50, columns 195 to 200, 5 pixels wide
bkgdMask(1:40,1:10) = true;
bkgdMask(30:50,195:200) = true;

% put variables in in struct, as getRawSignalsLineBkgdSub() expects
in.nMasks = 2;
in.nFrames = nFrames;
in.alSeries = fullImgSeriesStripes;
in.bkMask = bkgdMask;







