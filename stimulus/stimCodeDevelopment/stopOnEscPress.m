% stopOnEscPress.m
%
% Function to stop data acquisition from NIDAQ upon press of esc key 
% Uses KbCheck from Psychtoolbox, as there doesn't appear to be another way
% of doing this in Matlab
%

function stopOnEscPress(src, event)
    escCode = zeros(1,256);
    escCode(27) = 1;
    [~, ~, keyCode, ~] = KbCheck;
    if isequal(escCode,keyCode)
        disp('Stopping acquisition');
        plot(event.TimeStamps, event.Data);
        src.stop(); 
    end
end