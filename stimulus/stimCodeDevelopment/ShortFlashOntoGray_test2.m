% ShortFlashOntoGray_test2.m
% 
% The slightly more modular version of ShortFlashOntoGray_test1.m
% Note: displayContrastCycle is a nested function.
% Same code as ShortFlashOntoGray_FullFieldFlash.m
%
% last update: 08.29.16

classdef ShortFlashOntoGray_test2 < Stimulus 
    properties 
        Out
        ParamList = {'TotalDuration', 'FlashDuration', 'FlashContrast', ...
            'GrayDuration', 'GrayContrast'};
        TotalDuration % seconds
        FlashDuration  
        FlashContrast
        GrayDuration
        GrayContrast
    end 
    methods
        % Constructor
        function obj = ShortFlashOntoGray_test2(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi)
            % Fetch variables 
           
            totalDur = obj.TotalDuration{1};  % total duration of stimulus
            flashDur = cell2mat(obj.FlashDuration); % vector of different flash durations
            grayDur = obj.GrayDuration{1};
            grayFrames = round(grayDur/ifi);       
            flashContrast = cell2mat(obj.FlashContrast);% flash contrast values 
            gray = obj.GrayContrast{1}; % grey contrast value
            N = obj.nEpochs; % number of different flash-to-gray epochs
            % number of flash-to-grey sequences to play
            numSeq = ceil(totalDur/(mean(flashDur)+grayDur));
            
            % Create a random sequence flash-to-grey epochs
%             rng('shuffle'); % seeds the random number generator based on current time
            % Set random number generator to produce same sequence of
            % random numbers upon each all to displayStim()
            rng('default');
            obj.Out.RandSeed = rng; % save the random number generator seed
            randseq = ceil(N*rand(numSeq, 1));

            % --- Psychtoolbox code --- *
            numFrames = totalDur/ifi; % flip window at every frame
            startTime = Screen('Flip', window); 
            tic % clock for timing code efficiency
            currFrame = 1; % clock for keep track of current frame

            % Play the stimulus until total duration is reached or when
            % the user presses a key
            for i = 1:numSeq 
                % Draw a random number 
                k = randseq(i); 
                % Get the contrast value and duration of flash
                contrast = flashContrast(k);
                flashFrames = round(flashDur(k)/ifi);
                
                % For the duration of the flash, draw the uniform fullfield
                % contrast 
                currFrame = displayContrastCycle(contrast, ...
                    flashFrames, currFrame); 
                currFrame = displayContrastCycle(gray, grayFrames, ...
                    currFrame);
            end 
            
            Screen('CloseAll');
            
            toc % clock
            fprintf('current frame = %f', currFrame); % check that currFrame = totalFrames
            
            % Draws a fullfield contrast repeatedly at the frame rate   
            % for a given number of frames, flips the window
            function currFrame = displayContrastCycle(contrast, ...
                    totalFrames, currFrame)
                    for j = 1:totalFrames
                        currFrame = currFrame + 1;
                        % Draw the stimulus 
                        Screen('FillRect', window, contrast, obj.screenDim)
                        % Draw the input to the photodiode 
                        Screen('FillRect', window, contrast, obj.pdDim);
                         % Flip the window
                        Screen('Flip', window, startTime + currFrame*ifi);   
                        pause(0)
                        
                        if KbCheck
                            return;
                        end                    
                    end 
              end 
         
        end 
        
    end 
end 
