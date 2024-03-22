% WhiteNoise_1D_FullFieldFlash_playback.m
% 
% Same as WhiteNoise_1D_FullFieldFlash.m except generates the stimulus
% sequence prior to displaying the stimulus rather than generating stimulus
% values live. Here, PTB also takes a screenshot of each frame as it is
% displayed. The stimulus screenshots are saved as stimMovie.
%
% last update: 09.07.16

classdef WhiteNoise_1D_FullFieldFlash_playback < Stimulus 
    properties 
        Out
        ParamList = {'EpochDuration', 'EpochContrast',...
            'RepeatRNGSeed'};
        EpochDuration
        EpochContrast
        RepeatRNGSeed
    end 
    methods
        % Constructor
        function obj = WhiteNoise_1D_FullFieldFlash_playback(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)
            % Method A:
            % for each frame or each flip
            %   if done with previous epoch, draw a random contrast value 
            %   display the full field contrast 
            %   update frame counter for current contrast (countdown)
            % Method B:
            % for total number of stimulus updates (contrast values), 
            %   draw a random contrast value
            %   display the contrast for the given duration 
            
            % Here, I've implemented Method A 
            
            % Fetch stimlus paramters
            contrast = cell2mat(obj.EpochContrast); % vector of contrast values
            duration = obj.EpochDuration{1}; % duration of 1 contrast value
            nValues = length(contrast);  % total number of contrast values
            pdContrast = 0; % photodiode contrast value 
            
            % user-estimated total stimulus length
            estFrames = ceil(stimDuration/ifi); 
            buffer = ceil(180/ifi); % 3 min buffer
            % number of frames stimulus will actually run for
            totalFrames = estFrames + buffer;
%           % vector to save actual stimulus played
            rawStim = zeros(totalFrames, 1); 
            pdOut = zeros(totalFrames, 1); 
%             nUpdates = ceil(stimDuration/duration); % length of contrast sequence 
           
            % Set random number generator 
            if obj.RepeatRNGSeed{1} == 1 % use same random sequence for every trial you run the stimulus
                rng('default');
            else 
                rng('shuffle'); % use a different random number sequence for each trial
            end 

            % Generate a sequence of random contrast values whose length is
            % determined by the total duration of the stimulus and the
            % update rate. We draw from a uniform distribution because we
            % have such few contrast values 
        framesPerEpoch = int16(duration/ifi); % number of frames per contrast    
        nUpdates = ceil(totalFrames/framesPerEpoch);
        % random sequence of indices of stimulus epochs 
        stimInd = ceil(nValues*rand(nUpdates, 1));
        % reconstruct the stimulus using the saved random generator seed
        stimSeq = contrast(stimInd);
        frameStimVals = repelem(stimSeq, framesPerEpoch); % stim value at each frame
        
        % these should match
%         length(frameVal)
%         totalFrames
        
%             currFrame = 1; % initialize frame counter 
        % Arrays for saving outputs of Screen('Flip')
        missed = zeros(totalFrames, 1);
        vblTimes = zeros(totalFrames, 1);
        stimOnsetTimes = zeros(totalFrames, 1);
                
        % Matrix for saving PTB screenshots of the stimulus as it is
        % displayed
        totalFrames = 500;
        picSize = 50; % pixel dimensions of the screenshots
        imageGray = zeros(picSize, picSize, 100);
            
            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window); 
            
            % MAIN LOOP over each stimulus frame
            % Keep playing the stimulus until the user presses a key on the keyboard
            for currFrame = 1:totalFrames
                % photodiode switches btw W and B upon each epoch
                pdContrast = double(~pdContrast); 
                % set current frame's contrast value
                c = frameStimVals(currFrame); 
                % Draw the stimulus
                % contrasts are corrected for 6-bit RGB 
                Screen('FillRect', window, c*obj.RGBscale, obj.screenDim); 
                % Draw onto the photodiode
                Screen('FillRect', window, pdContrast*obj.RGBscale, obj.pdDim);
                
                % save the contrast of most recent frame
                rawStim(currFrame) = c; 
                pdOut(currFrame) = pdContrast;
                
                % Flip window at the frame rate
                [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed]...
                    = Screen('Flip', window, startTime + currFrame*ifi);  
                pause(0)
                
                imageRGB = Screen('GetImage', window, [200 200 200+picSize 200+picSize]);
                imageGray(:,:,currFrame) = rgb2gray(imageRGB);
                
                 % Test code
                missed(currFrame) = Missed;
                vblTimes(currFrame) = VBLTimestamp;
                stimOnsetTimes(currFrame) = StimulusOnsetTime;
                
                % If user hits keyboard, exit the stimulus loop 
                if KbCheck
                    currFrame
                    break;
                end    
               
%                 currFrame = currFrame + 1; % update frame counter
            end
            
            Screen('CloseAll');
            
            % Save metadata specific to this class
            r = rng;
            obj.Out.rndSeed = r.Seed; % seed for random number generator  
            obj.Out.rawStim = rawStim(1:currFrame); % save stimulus contrast values to Out 
            obj.Out.frameStimVals = frameStimVals(1:currFrame);
            
            % Code for timing tests
            obj.Out.pdOut = pdOut(1:currFrame); % Test code
            obj.Out.Missed = missed(1:currFrame);
            obj.Out.VBLTimestamp = vblTimes(1:currFrame);
            obj.Out.StimulusOnsetTime = stimOnsetTimes(1:currFrame);
            obj.Out.stimMovie = imageGray;
        end 
        
    end 
end 
