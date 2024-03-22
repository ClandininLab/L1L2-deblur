% computeFramePeaks.m
%
% Function that takes in average impulse response trace and returns the
%  index where the first phase peaks and the index where the second phase
%  peaks
%
% INPUT:
%   avgTrace - the impulse response trace
%   whichContrast - whether it's a light or dark impulse response (1 for 
%       dark and 2 for light)
%
% OUTPUT:
%   framePeak1 - index where first phase peaks
%   framePeak2 - index where second phase peaks
%
% Updated: HHY 10/12/17
%

function [framePeak1, framePeak2] = computeFramePeaks(avgTrace, ...
    whichContrast)
    
    if (whichContrast == 1) % dark
        % dark = depol = negative = min
        [~, framePeak1] = min(avgTrace); 
        % second peak must be after first
        [~, peak2] = max(avgTrace(framePeak1:end));
        framePeak2 = peak2 + framePeak1 - 1;
    elseif (whichContrast == 2) % light
        % opposite for light trace
        [~, framePeak1] = max(avgTrace);
        [~, peak2] = min(avgTrace(framePeak1:end));
        framePeak2 = peak2 + framePeak1 - 1;
    else
        framePeak1 = 0;
        framePeak2 = 0;
    end
end