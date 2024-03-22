% correctFFFBinNaN.m
%
% Function to correct any NaNs in binning of full field flash traces.
%  Replace with mean of adjacent 2 bins or if adjacent bins also contain
%  NaNs, replace with median of trace. Std Err set to 0 in this case
%
% INPUTS:
%     corTrace - trace to be corrected for NaNs
%     wrapTrace - trace to use for edge effects end preceeds beginning of
%        corTrace and beginning follows end of corTrace
%     corStdErr - stdErr of trace to be corrected
%
% OUTPUTS:
%     meanTrace - corTrace corrected for NaNs, if no NaNs, equals corTrace
%     stdErrs - stdErr trace, 0 for NaNs, if no NaNs, equals corStdErr
%

function [meanTrace, stdErrs] = correctFFFBinNaN(corTrace, wrapTrace,...
    corStdErr)

    corNaN = find(isnan(corTrace));
    
    meanTrace = corTrace;
    stdErrs = corStdErr;
    
    for i = 1:length(corNaN)
        % deal with wrap around
        if (corNaN(i)==1)
            leftVal = wrapTrace(end);
            rightVal = corTrace(corNaN(i)+1);
        elseif (corNaN(i)==length(corTrace))
            rightVal = wrapTrace(1);
            leftVal = corTrace(corNaN(i)-1);
        else
            leftVal = corTrace(corNaN(i)-1);
            rightVal = corTrace(corNaN(i)+1);
        end

        % mean of adjacent bins
        if (~isnan(leftVal) && ~isnan(rightVal))
            meanTrace(corNaN(i)) = mean([leftVal rightVal]);
        else % median of whole trace (non-NaN values
            meanTrace(corNaN(i)) = median(corTrace(~isnan(corTrace)));
        end
        
        stdErrs(corNaN(i)) = 0;
    end
end