% WhiteNoise_1D_FullFieldFlash_test.m
% 
% Same as WhiteNoise_1D_FullFieldFlash.m except PTB takes a screenshot of each
% stimulus frame displayed during stimulus presentation. 
%
% last update: 09.07.16

classdef WhiteNoise_1D_FullFieldFlash_test < Stimulus 
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
        function obj = WhiteNoise_1D_FullFieldFlash_test(txtFile)
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
            epochFrameCount = 0; % frame counter to track frames left in an epoch
            
            % user-estimated total stimulus length
            estFrames = ceil(stimDuration/ifi); 
            buffer = ceil(180/ifi); % 3 min buffer
            % number of frames stimulus will actually run for
            totalFrames = estFrames + buffer;
%           % vector to save actual stimulus played
            rawStim = zeros(totalFrames, 1); 
            pdOut = zeros(totalFrames, 1); 
            nUpdates = 0; % counts number of stimulus contrast updates
           
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
%             randseq = ceil(nValues*rand(nUpdates+1, 1)); 
%             i = 1;
        
%             currFrame = 1; % initialize frame counter 

            % Arrays for saving outputs of Screen('Flip')
            missed = zeros(totalFrames, 1);
            vblTimes = zeros(totalFrames, 1);
            stimOnsetTimes = zeros(totalFrames, 1);
            
            % Matrix for saving PTB screenshots of the stimulus as it is
            % displayed
            totalFrames = 1000;
            picSize = 50; % pixel dimensions of the screenshots
            imageGray = zeros(picSize, picSize, totalFrames);
            

            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window); 
            
            % MAIN LOOP over each stimulus frame
            % Keep playing the stimulus until the user presses a key on the keyboard
            for currFrame = 1:totalFrames 
                % If the previous epoch ended, pick a new contrast value to 
                % start a new epoch
                if  epochFrameCount == 0 
%                   k = randseq(i); % use if random sequence was generated prior to for loop
%                   i = i + 1;
                    k = ceil(nValues*rand(1)); % use if drawing random number upon each frame
                    c = contrast(k);
                    pdContrast = double(~pdContrast); % photodiode switches btw W and B upon each epoch
                    epochFrameCount = round(duration/ifi); % load the counter 
                    nUpdates = nUpdates + 1; % update stimulus epoch counter
                end 
                
                % Draw the stimulus
                % contrasts are corrected for 6-bit RGB 
                Screen('FillRect', window, c*obj.RGBscale, obj.screenDim); 
                % Draw onto the photodiode
                Screen('FillRect', window, pdContrast*obj.RGBscale, obj.pdDim);
                
%                 % testing code: used this when photodiode was positioned in front
%                 % of the screen 
%                 if (c == 1 || c == 0)
%                     Screen('FillRect', window, c, obj.screenDim); 
%                 else 
%                     Screen('FillRect', window, [1 0 1], obj.screenDim);
%                 end 

%                 % test: 
%                 if (c < 1 && c > 0)
%                     Screen('FillRect', window, [1 0 1], obj.pdDim);
%                 else
%                     Screen('FillRect', window, pdContrast, obj.pdDim);
%                 end 
                
                % Flip window at the frame rate
                [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed]...
                    = Screen('Flip', window, startTime + currFrame*ifi);  
                pause(0)
                
                % Take a screenshot of the stimulus frame displayed
                imageRGB = Screen('GetImage', window, [200 200 200+picSize 200+picSize]);
                imageGray(:,:,currFrame) = rgb2gray(imageRGB);
                
                % save the contrast of most recent frame
                rawStim(currFrame) = c; 
                pdOut(currFrame) = pdContrast; % Test code
                missed(currFrame) = Missed;
                vblTimes(currFrame) = VBLTimestamp;
                stimOnsetTimes(currFrame) = StimulusOnsetTime;
                
                % update frame counter for this stimulus epoch
                epochFrameCount = epochFrameCount - 1;
                
                % If user hits keyboard, exit the stimulus loop 
                if KbCheck
                    break;
                end    
               
%                 currFrame = currFrame + 1; % update frame counter
            end
            
            Screen('CloseAll');
            
            % Save metadata specific to this class
            r = rng;
            obj.Out.rndSeed = r.Seed; % seed for random number generator  
            obj.Out.rawStim = rawStim(1:currFrame); % save stimulus contrast values to Out 
            obj.Out.nUpdates = nUpdates;
            
            % Code for timing tests
            obj.Out.pdOut = pdOut(1:currFrame); % Test code
            obj.Out.Missed = missed(1:currFrame);
            obj.Out.VBLTimestamp = vblTimes(1:currFrame);
            obj.Out.StimulusOnsetTime = stimOnsetTimes(1:currFrame);
            obj.Out.stimMovie = imageGray;
        end 
        
    end 
end 
