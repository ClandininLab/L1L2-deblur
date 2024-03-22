% generate_sineMonophasicFilter.m
%
% Function to generate filter similar to those observed for bipolar cells.
% First portion of sine wave
%
% INPUT:
%   startInd - index where filter starts (zero before this)
%   peak1Val - value of filter's first peak
%   filterLength - length of filter, not including zeros
%   totFiltLength - total length of filter in sample points, including
%       zeros
%
% OUTPUT:
%   filter - row vector of length filtLength representing filter
%
% Updated: 2/10/18 - Helen Yang
%

function filter = generate_sineMonophasicFilter(startInd,  ...
    peak1Val, filterLength, totFiltLength)

    filter = zeros(1,totFiltLength);
    
    % conversion factor between indicies and degrees; 180 as half sine wave
    d2i = filterLength/180;
    
    % indicies in degrees
    filtPoints = 0:(1/d2i):floor(filterLength/d2i);
    
    % get filter
    filter(startInd:(startInd + length(filtPoints)-1)) = ...
        peak1Val * sind(filtPoints);
    
end