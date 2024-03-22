% FullFieldFlash_NEpochs.m
%
% *****DOES NOT WORK YET!!!!!!!!!***
% 
% An "epoch" is defined as a contrast value and a duration (seconds).
% User may specify any number of epochs of any duration.
% Iterates through each epoch in the .txt file sequentially for amount of
% time specified by user.
%
% This code was derived from FullFieldFlash_BinaryContrast.m and
% SearchStimulusFlash.m
% 
% Parameters:  
%       Duration - duration of each contrast epoch in seconds
%       Contrast - contrast values for each epoch 
%
% Input: (see '2sFFF_template.txt' as an example)
%       txtFile - .txt file specifying stimulus parameters and stimulus
%           class name
%
% Output (stim.Out):
%       rawStim - vector containing the stimulus contrast value played at 
%           each frame
%
% last update: 09.09.16

classdef FullFieldFlash_NEpochs < Stimulus 
    properties 
        Out
        ParamList = {'Duration', 'Contrast'};
        Duration
        Contrast
    end 
    methods
        % Constructor initializes an object of this stimulus class using
        % the .txt file chosen by the user
        function obj = FullFieldFlash_NEpochs(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)  
            % estimate of total number of stimulus frames based on duration
            % specified by user
            estFrames = ceil(stimDuration/ifi);
            buffer = ceil(60/ifi); % 1min buffer
            totalFrames = estFrames + buffer;
%           % vector to save actual stimulus played
            rawStim = zeros(totalFrames, 1); 
            
            epochCount = 1;
            
            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window);
            
            % keep playing the stimulus until the user presses a key on the keyboard
            for frame = 1:totalFrames 
                currEpoch = mod(epochCount-1,obj.nEpochs)+1;
                epochCount = ceil((frame*ifi)/obj.Duration{currEpoch});
             
%                 if (mod(currEpoch, 2)) % light epochs are when currEpoch is odd
%                     contrast = GrayIndex(window, lightContrast);
%                 else  % dark epochs are even
%                     contrast = GrayIndex(window, darkContrast);
%                 end 
                contrast = obj.Contrast{currEpoch};
                
                % Draw the stimulus 
                % contrasts are corrected for 6-bit RGB coding
                Screen('FillRect', window, contrast*obj.RGBscale, obj.screenDim); 
                % Draw the input to the photodiode 
                Screen('FillRect', window, contrast*obj.RGBscale, obj.pdDim);
                
                 % flip entire window at the frame rate
                Screen('Flip', window, startTime + frame*ifi);  
                
                pause(0) % super important - without it, NIDAQ won't scan
                rawStim(frame) = contrast; 
                
                % If user hits keyboard, exit the stimulus loop 
                if KbCheck
                    break;
                end    
            end 
            
            Screen('CloseAll');
            
            rawStim = rawStim(1:frame);
            obj.Out.rawStim = rawStim;
        end 
        
    end 
end 
