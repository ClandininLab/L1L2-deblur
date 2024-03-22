% FullFieldFlash_BinaryContrast_RGB.m
% 
% Version of FullFieldFlash_BinaryContrast that works with lightcrafter
%  pattern setup where 3 patterns are different (the same one as in 
%  NaturalisticStimulus_1D_FullField). Note that as currently written, the
%  flash durations have the same temporal resolution as in 
%  FullFieldFlash_BinaryContrast and not 3X the resolution, though this
%  could be changed if desired.
% Subclass of Stimulus that alternates between 2 epochs, each defined by a 
%  contrast value and a duration. Presents that contrast over full screen 
%  for that duration.
%  Code may be revised to include N epochs. 
%
% This code was derived from FullFieldFlash_BinaryContrast.
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
%   8/9/19 - HHY - 
% last update: 8/9/19 HHY
% 
classdef FullFieldFlash_BinaryContrast_RGB < Stimulus 
    properties 
        Out
        ParamList = {'Duration', 'Contrast'};
        Duration
        Contrast
    end 
    methods
        % Constructor initializes an object of this stimulus class using
        % the .txt file chosen by the user
        function obj = FullFieldFlash_BinaryContrast_RGB(txtFile)
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
                
                % convert contrasts for RGB encoding
                if (mod(epochCount, 2)) % light epochs are when currEpoch is odd
                    contrast = convertInt2RGB(...
                        repmat(GrayIndex(window, lightContrast), 1, 3));
                    pdContrast = convertInt2RGB(ones(1, 3));
                else  % dark epochs are even
                    contrast = convertInt2RGB(...
                        repmat(GrayIndex(window, darkContrast), 1, 3));
                    pdContrast = convertInt2RGB(zeros(1, 3));
                end 
                
                % Draw the stimulus 
                Screen('FillRect', window, contrast, obj.screenDim); 
                % Draw the input to the photodiode 
                Screen('FillRect', window, pdContrast, obj.pdDim);
                
                 % flip entire window at the frame rate
                [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed]...
                    = Screen('Flip', window, startTime + frame*ifi);  
                
                pause(0) % super important - without it, NIDAQ won't scan
                % save contrast of frame, converted back to intensity from
                %  RGB; take first since all same
                thisFrameContrast = convertRGB2Int(contrast);
                rawStim(frame) = thisFrameContrast(1); 
                
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
