% computeFrameZero2.m
%
% Function that takes in average impulse response trace and returns the
%  zero crossing between the first and second phases
%
% INPUT
%   avgTrace - the average impulse response trace
%   framePeak - index of peak of 1st phase 
%   whichContrast - whether trace is light or dark impulse response; 1 for
%       dark and 2 for light
%
% OUPUT
%   frameZero - index of frame where zero crossing happens (first frame on
%       other side of zero)
%
% Updated: 10/12/17

function frameZero = computeFrameZero2(avgTrace, framePeak, whichContrast)

    % dark - zero crossing is negative to positive 
    if (whichContrast == 1)
        % zero crossing frame is first frame on other side of zero after
        %  peak of 1st phase; for dark, first positive value
        frameZero = find(avgTrace(framePeak:end) >= 0, 1, 'first') + ...
            framePeak - 1; % add framePeak-1 to account for index diff
    elseif (whichContrast == 2) % light
        % first negative value
        frameZero = find(avgTrace(framePeak:end) <= 0, 1, 'first') + ...
            framePeak - 1;       
    else
        frameZero = 0;
    end

end