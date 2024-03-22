% ShortFlashOntoGray_test3.m
% 
%  Most modular/efficient version of ShortFlashOntoGray_ so far. 
%
% last update: 09.06.16

classdef ShortFlashOntoGray_test3 < Stimulus 
    properties 
        Out
        ParamList = {'FlashDuration', 'FlashContrast', ...
            'GrayDuration', 'GrayContrast', 'RepeatRNGSeed'};
        FlashDuration  
        FlashContrast
        GrayDuration
        GrayContrast
        RepeatRNGSeed
    end 
    methods
        % Constructor
        function obj = ShortFlashOntoGray_test3(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)
            
            % Fetch parameters 
            totalDur = stimDuration;  % total duration of stimulus in seconds
            flashDur = cell2mat(obj.FlashDuration); % vector of different flash durations
            grayDur = obj.GrayDuration{1}; % duration of gray epoch in seconds 
            totalFrames = round(totalDur/ifi); 
            flashFrames = round(flashDur./ifi);  % vector of #frames for flash types
            grayFrames = round(grayDur/ifi);       
            flashContrast = cell2mat(obj.FlashContrast); % flash contrast values 
            gray = obj.GrayContrast{1}; % grey contrast value

            gray = obj.RGBscale*gray;
            nFlashType = length(flashContrast); % number of different flash types
            
            % Set random number generator 
            if obj.RepeatRNGSeed{1} == 1 
                % Set random number generator to produce same sequence of
                % random numbers upon each all to displayStim()
                rng('default');
            else 
                rng('shuffle'); % use a different random number sequence for each trial
            end 
            
            % Vector to save stimulus contrast value at each frame
            rawStim = zeros(totalFrames, 1);

            % Counters to keep track of current position in flash-to-gray sequence 
            epochFrameCount = 0;
            isGray = 0;

            startTime = Screen('Flip', window); 

            for currFrame = 1:totalFrames
                % If previous gray epoch ended, randomly draw a new flash epoch 
                if epochFrameCount == 0 && isGray == 1 
                    k = ceil(nFlashType*rand(1));
                    contrast = flashContrast(k);
                    epochFrameCount = flashFrames(k); % run new contrast value for this many frames
                    isGray = 0; % switch to flash epoch 

                % If previous flash period ended, start grey epoch
                elseif epochFrameCount == 0 && isGray == 0
                    contrast = gray;
                    epochFrameCount = grayFrames; 
                    isGray = 1; 
                end 
                
%                 % Test the timing of the epochs (comment out later)
%                 if ~(contrast == 1 || contrast == 0) % grey
%                        Screen('FillRect', window, 0.5 , obj.screenDim);
% %                      Screen('FillRect', window, [1 0 1], obj.screenDim); % works if RGB pattern sequence works
%                 else
%                      Screen('FillRect', window, contrast, obj.screenDim); % stimulus screen
%                 end 

                % set contrast for photodiode  
                if (contrast == gray) % use 0.5 for gray for better visibility 
                    pdContrast = 0.5;
                else
                    pdContrast = contrast;
                end
                    
                % Draw full field flash with the current contrast value 
                Screen('FillRect', window, contrast, obj.screenDim); % stimulus screen
                Screen('FillRect', window, pdContrast, obj.pdDim); % photodiode input
                % Flip the window at the framerate 
                Screen('Flip', window, startTime + currFrame*ifi);   
                pause(0)
                
                % save the contrast value at each frame
                rawStim(currFrame) = contrast;
                % update frame counter for this stimulus epoch
                epochFrameCount = epochFrameCount - 1;
                
                % Exit this function if user hits keyboard. Does not save
                % stimulus data!
                if KbCheck
                    break;
                end       
        
            end 

            Screen('CloseAll');
            
            % save metadata specific to this class
            r = rng;
            obj.Out.rndSeed = r.Seed; % seed for random number generator 
            obj.Out.rawStim = rawStim; % save stimulus contrast values to Out
            
        end 
    end 
end 