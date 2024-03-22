function [time] = teststim_FFF(stim, stimDuration, screenNumber, window, ifi)
% Full Field Flash
% Use playstim_v2.m to test me.

% contrasts
black = stim.contrast{1}; % right now, it's just 1 or 0, we will make this more rigorous
white = stim.contrast{2};

% seconds per stimulus epoch, assuming that light and dark duration are the
% same length
epochDur = stim.duration{1}; 
% waitframes = round(epochDur / ifi);

numFlips = ceil(stimDuration/ifi);

startTime = Screen('Flip', window);
% startTime = GetSecs();
% i=1;
for i=1:numFlips
    currEpoch = floor((i*ifi)/epochDur);
    if (mod(currEpoch,2))
        Screen('FillRect', window, white);
        Screen('Flip', window, startTime + i*ifi);
    else 
        Screen('FillRect', window, black);
        Screen('Flip', window, startTime + i*ifi);
    end    
    pause(0)
    
    % keyboard press also stops stimulus
    if KbCheck 
        sca;
        break;
    end
end 

time = GetSecs();

% % Clear the screen.
% sca;
end 
