% generate_peakyBiphasicFilter.m
%
% Function to generate filter similar to those observed for L1 and L2 at
%  baseline. Linear from start to peak1, linear from peak1 to peak2, single
%  exponential decay. Filter is OFF signed. Units are sample points, not
%  time.
%
% INPUT:
%   startInd - index where filter starts (zero before this)
%   peak1Ind - index where filter first peaks
%   peak2Ind - index where filter's 2nd phase peaks
%   peak1Val - value of filter's first peak
%   peak2Val - value of filter's second peak
%   lambda - decay constant of single exponential decay from 2nd peak
%   filtLength - length of filter in sample points
%
% OUTPUT:
%   filter - row vector of length filtLength representing filter
%
% Updated: 2/9/18 - Helen Yang
%

function filter = generate_peakyBiphasicFilter(startInd, peak1Ind, ...
    peak2Ind, peak1Val, peak2Val, lambda, filtLength)

    filter = zeros(1,filtLength);
    
    % decay of second phase
    for i = peak2Ind:length(filter)
        ind = i - peak2Ind;
        filter(i) = peak2Val * exp(-1/lambda * ind);
    end

    % 1st phase to peak
    slope1 = peak1Val / (peak1Ind - startInd);
    for j = startInd:peak1Ind
        ind = j-startInd;
        filter(j) = slope1 * ind;
    end

    % 1st phase peak to 2nd phase peak
    slope2 = (peak2Val-peak1Val)/(peak2Ind-peak1Ind);
    for k = peak1Ind:peak2Ind
        ind = k - peak1Ind;
        filter(k) = slope2 * ind + peak1Val;
    end 
end