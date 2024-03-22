% FullFieldFlashSubclass is a subclass of the Stimulus superclass.
% last update: 08.19.16

classdef FullFieldFlashSubclass < Stimulus 
    properties 
        Out
        ParamList = {'Duration', 'Contrast'};
        Duration
        Contrast
        
    end 
    methods
        % Constructor
        function obj = FullFieldFlashSubclass(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi)
            % Contrast values for each epoch. Currently either 1 or 0,
            % we'll make this more rigorous in coming versions
            darkContrast = obj.Contrast{1};
            lightContrast = obj.Contrast{2};
            
            % Seconds per stimulus epoch, assuming that light and dark
            % durations are of equal length
            epochDur = obj.Duration{1};
            
            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window);
            
            frame = 0;
            % keep playing the stimulus until the user presses a key on the keyboard
            while ~KbCheck 
                currEpoch = floor((frame*ifi)/epochDur);
                % light epoch
                if (mod(currEpoch, 2))
                    Screen('FillRect', window, lightContrast); % draw the stimulus
                    Screen('Flip', window, startTime + frame*ifi); % flip window at the frame rate
                % dark epoch
                else
                    Screen('FillRect', window, darkContrast);
                    Screen('Flip', window, startTime + frame*ifi);  
                end 
                pause(0) % super important - without it, NIDAQ won't scan
                frame = frame + 1;
            end 
         
        end 
        
    end 
end 
