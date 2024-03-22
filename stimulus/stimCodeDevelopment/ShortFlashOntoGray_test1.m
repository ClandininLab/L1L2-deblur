% ShortFlashOntoGray_FullFieldFlash.m
% 
% displayStim() in this version is the "long" implementation 
% where modularity is not attempted 
%
% last update: 08.29.16

classdef ShortFlashOntoGray_test1 < Stimulus 
    properties 
        Out
        ParamList = {'TotalDuration', 'FlashDuration', 'FlashContrast', ...
            'GrayDuration', 'GrayContrast'};
        TotalDuration
        FlashDuration
        FlashContrast
        GrayDuration
        GrayContrast
    end 
    methods
        % Constructor
        function obj = ShortFlashOntoGray_test1(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi)          
            % Grab variables 
            % total duration of stimulus
            totalDur = obj.TotalDuration{1}; 
            % vector of different flash durations
            flashDur = cell2mat(obj.FlashDuration);
            grayDur = obj.GrayDuration{1};
            grayFrames = round(grayDur/ifi);
            % contrast values 
            flashContrast = cell2mat(obj.FlashContrast);
            gray = obj.GrayContrast{1};
            % number of different types of flash-to-gray epochs
            N = obj.nEpochs; 
            % number of flash-to-grey sequences to play
            numSeq = ceil(totalDur/(mean(flashDur)+grayDur));
            
            % Create a random sequence flash-to-grey epochs
            % Set random number generator to produce same sequence of
            % random numbers upon each all to displayStim()
            rng('default');
            obj.Out.RandSeed = rng; % save the random number generator seed
            randseq = ceil(N*rand(numSeq, 1));

            % --- Psychtoolbox code --- *
            numFrames = totalDur/ifi % total number of frames 
            startTime = Screen('Flip', window);      
            tic % clock
            currFrame = 1; % clock to keep track of current frame 
            
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
                for j = 1:flashFrames
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
                
                contrast = gray;
                % Play the gray interweave 
                for j = 1:grayFrames
                    currFrame = currFrame + 1;
                    % Draw the stimulus 
                    Screen('FillRect', window, contrast, obj.screenDim)
                    % Draw the input to the photodiode 
                    Screen('FillRect', window, contrast, obj.pdDim);
                     % Flip the window
                    Screen('Flip', window, startTime+ currFrame*ifi);  
                    pause(0)              
                    if KbCheck
                        return;
                    end
                end 


            end 
            
            toc % clock
            
            sca;
            
            fprintf('current frame = %f', currFrame); % check that currFrame = totalFrames
         
        end 
        
    end 
end 
