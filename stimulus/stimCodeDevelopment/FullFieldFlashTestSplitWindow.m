% FullFieldFlashTestSplitWindow.m
% 
% Same as FullFieldFlashSubclass except displays the stimulus to a small
% region in the screen window to the photodiode in addition to the actual
% stimulus.
%
% last update: 08.29.16

classdef FullFieldFlashTestSplitWindow < Stimulus 
    properties 
        Out
        ParamList = {'Duration', 'Contrast'};
        Duration
        Contrast
    end 
    methods
        % Constructor
        function obj = FullFieldFlashTestSplitWindow(txtFile)
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
            
%             % vector to save actual stimulus played
            rawStim = zeros(30000, 1); % saves data up to 5min 
            
            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window);
            
            frame = 1;
            % keep playing the stimulus until the user presses a key on the keyboard
            while ~KbCheck 
                currEpoch = ceil((frame*ifi)/epochDur);
             
                if (mod(currEpoch, 2)) % light epochs are when currEpoch is odd
                    contrast = GrayIndex(window, lightContrast);
                else  % dark epochs are even
                    contrast = GrayIndex(window, darkContrast);
                end 
                
                % Draw the stimulus 
                Screen('FillRect', window, contrast, obj.screenDim); 
                % Draw the input to the photodiode 
                Screen('FillRect', window, contrast, obj.pdDim);
                
                 % flip entire window at the frame rate
                Screen('Flip', window, startTime + frame*ifi);  
                
                pause(0) % super important - without it, NIDAQ won't scan
                rawStim(frame) = contrast; 
                frame = frame + 1;
            end 
            
            Screen('CloseAll');
            
            rawStim = rawStim(1:frame);
            obj.Out.rawStim = rawStim;
        end 
        
    end 
end 
