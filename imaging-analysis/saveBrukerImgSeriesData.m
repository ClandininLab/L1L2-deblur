% saveBrukerImgSeriesData.m
%
% Saves time series data into time series folder as file imDat.mat
%
% INPUT:
%   seriesName - name of time series
%   alignedSeries - image sequence, aligned time series
%   unalignedSeries - image sequence, original, unaligned time series
%   meanImageAligned - mean over time of the whole image sequence
%   xDim - size of imaging frame, x dimension
%   yDim - size of imaging frame, y dimension
%   frameTimes - frame times of images
%   tSeriesPath - full path to TSeries folder
%
% OUTPUT:
%   no output, saves time series data as imDat.mat file into TSeries folder
%    as imDat struct
%

function saveBrukerImgSeriesData(seriesName, alignedSeries, ...
    unalignedSeries, meanImageAligned, xDim, yDim, frameTimes, tSeriesPath)

    % image data for this time series
    imDat = struct; 
    imDat.fileloc = tSeriesPath; % file location
    imDat.flyName = ''; % maybe update this
    imDat.seriesName = seriesName; % series name
    imDat.alSeries = alignedSeries; % aligned time series
    imDat.unSeries = unalignedSeries; % unaligned series
    imDat.avgIm = meanImageAligned; % average image
    imDat.nFrames = size(alignedSeries, 3);
    imDat.xPixels = xDim; % pixels in x dimension of FOV
    imDat.yPixels = yDim; % pixels in y dimension of FOV
    imDat.imFrameStartTimes = frameTimes; % frame times, in sec
    imDat.imIFI = mean(diff(frameTimes)); % ifi, in sec
    
    imDatFilepath = [tSeriesPath filesep 'imDat.mat'];

    save(imDatFilepath, '-struct', 'imDat');
end