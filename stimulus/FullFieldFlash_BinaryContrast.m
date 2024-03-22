% FullFieldFlash_BinaryContrast.m
% 
% Subclass of Stimulus that alternates between 2 epochs, each defined by a 
%  contrast value and a duration. Presents that contrast over full screen 
%  for that duration.
%  Code may be revised to include N epochs. 
%
% This code was derived from FullFieldFlashTestSplitWindow.
% 
% Parameters:  **NOTE: max number of epochs is 2 **
% 	Duration - duration of each contrast epoch in seconds
% 	Contrast - contrast values for each epoch 
%
% Constructor: obj = FullFieldFlash_BinaryContrast(txtFile) 
% 	txtFile - .txt file specifying stimulus parameters and stimulus
%   	class name (see '2sFFF_template.txt' as an example)
%
% Output (adds to stim.Out):
% 	rawStim - vector containing the stimulus contrast value played at 
%       each frame
%   stimIFI - IFI of stimulus presentation, inherited from psychtoolbox
%   Missed - psychtoolbox readout of missed frames (for debugging)
%   VBLTimestamp - psychtoolbox VBL time stamps for each frame (for
%       debugging)
%   StimulusOnsetTime - psychtoolbox output of StimulusOnsetTime (for
%       debugging)
%
% Updates:
%   3/7/17 - added definition of reconstructStim method
%   3/8/17 - added definition of selectDFFbaseline method
%   3/27/19 - corrected so photodiode signal is still full black/white even
%       if contrast of stimulus on screen is not
% last update: 3/27/19 HHY
% 
classdef FullFieldFlash_BinaryContrast < Stimulus 
    properties 
        Out
        ParamList = {'Duration', 'Contrast'};
        Duration
        Contrast
    end 
    methods
        % Constructor initializes an object of this stimulus class using
        % the .txt file chosen by the user
        function obj = FullFieldFlash_BinaryContrast(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)  
            
            % If more than 2 contrast values, print error message and exit
            % function. 
            if length(obj.Contrast) > 2
                display('Too many epochs in the stimulus. Please revise .txt file and try again.');
                return;
            end 
            
            % Contrast values for each epoch 
            darkContrast = (obj.Contrast{1});
            lightContrast = (obj.Contrast{2});
            
            % Seconds per stimulus epoch, assuming that light and dark
            % durations are of equal length
            epochDur = obj.Duration{1};
            % estimate of total number of stimulus frames based on duration
            % specified by user
            estFrames = ceil(stimDuration/ifi);
            buffer = ceil(180/ifi); % 3min buffer
            totalFrames = estFrames + buffer;
%           % vector to save actual stimulus played
            rawStim = zeros(totalFrames, 1); 
            
            % Arrays for saving outputs of Screen('Flip')
            missed = zeros(totalFrames, 1);
            vblTimes = zeros(totalFrames, 1);
            stimOnsetTimes = zeros(totalFrames, 1);

            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window);
            
            % keep playing the stimulus until the user presses a key on the keyboard
            for frame = 1:totalFrames 
                epochCount = ceil((frame*ifi)/epochDur);
                
                if (mod(epochCount, 2)) % light epochs are when currEpoch is odd
                    contrast = GrayIndex(window, lightContrast);
                    pdContrast = 1;
                else  % dark epochs are even
                    contrast = GrayIndex(window, darkContrast);
                    pdContrast = 0;
                end 
                
                % Draw the stimulus 
                % contrasts are corrected for 6-bit RGB coding
                Screen('FillRect', window, contrast*obj.RGBscale, obj.screenDim); 
                % Draw the input to the photodiode 
                Screen('FillRect', window, pdContrast*obj.RGBscale, obj.pdDim);
                
                 % flip entire window at the frame rate
                [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed]...
                    = Screen('Flip', window, startTime + frame*ifi);  
                
                pause(0) % super important - without it, NIDAQ won't scan
                rawStim(frame) = contrast; 
                
                % save PTB timing information
                missed(frame) = Missed;
                vblTimes(frame) = VBLTimestamp;
                stimOnsetTimes(frame) = StimulusOnsetTime;
                
                % If user hits keyboard, exit the stimulus loop 
                if KbCheck
                    break;
                end    
            end 
            
            Screen('CloseAll');
            
            rawStim = rawStim(1:frame);
            obj.Out.rawStim = rawStim;
            obj.Out.stimIFI = ifi;
            obj.Out.Missed = missed(1:frame);
            obj.Out.VBLTimestamp = vblTimes(1:frame);
            obj.Out.StimulusOnsetTime = stimOnsetTimes(1:frame);
        end 
        
        % For reconstructing the stimulus
        function [rcStim, rcStimInd] = reconstructStim(obj, ...
        	stimEpochStartTimes, lightStartTimes, darkStartTimes, order)
            
            display('Processing photodiode data for FullFieldFlash_BinaryContrast ...');

            % reconstruct the stimulus
            % preallocate
            rcStimInd = zeros(length(stimEpochStartTimes), 1);
            
            % determine whether timestamps start on dark or light epoch
            if min(darkStartTimes) < min(lightStartTimes)
                for c = 1:length(stimEpochStartTimes)
                    rcStimInd(c) = mod(c-1, obj.nEpochs)+1;
                end
            else
                for c = 1:length(stimEpochStartTimes)
                    rcStimInd(c) = mod(c, obj.nEpochs)+1;
                end 
            end
            
            contrasts = cell2mat(obj.Contrast);
            rcStim = contrasts(rcStimInd)';        
        end
        
        % selects imaging frames to be used as F0 in computing dF/F during
        %  analysis; here, all of them
        function [baselineFrameTimes, baselineSignals] = ...
                selectDFFbaseline(obj,imgFrameTimes, bksSignal, ...
                lightStartTimes, darkStartTimes)    
        	baselineFrameTimes = imgFrameTimes;
        	baselineSignals = bksSignal;
        end    
    end
end 
