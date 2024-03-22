% computeFrameZero1.m
%
% Function that takes in average impulse response trace and returns the
%  first zero crossing (or more precisely, the inflection point at the
%  start of the trace)
%
% INPUT
%   avgTrace - the average impulse response trace
%   framePeak - index of peak of 1st phase 
%   whichContrast - whether trace is light or dark impulse response; 1 for
%       dark and 2 for light
%
% OUPUT
%   frameZero - index of frame where zero crossing happens
%
% Updated: 10/12/17

function frameZero = computeFrameZero1(avgTrace, framePeak, whichContrast)
    % assume signal at 0 before stimulus started
    avgTrace = [0 avgTrace];
    % get difference
    diffTrace = diff(avgTrace);
    
    frameZero = 2; % default for zero crossing is frame 2
    for i = 1:framePeak % has to happen before peak of 1st phase
        if (whichContrast == 1) % dark
            % difference between frame and previous one is greater than 10%
            %  plus previous change (less than because dark negative)
            if (diffTrace(i+1) < (diffTrace(i) + 0.1*diffTrace(i+1)))
                frameZero = i;
                break;
            end
        elseif (whichContrast == 2) % light
            if (diffTrace(i+1) > (diffTrace(i) + 0.1*diffTrace(i+1)))
                frameZero = i;
                break;
            end            
        else
            frameZero = 0;
        end
    end
end