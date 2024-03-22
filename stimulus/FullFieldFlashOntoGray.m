% FullFieldFlashOntoGray.m
% 
% Subclass of Stimulus that presents a flash followed by gray, both over
%  the full extent of the screen.
% A 'sequence' is defined as a full field contrast ('Flash' period) 
%  followed by a gray interleave. This stimulus class accepts any number of 
%  Flash types, each defined by a contrast value between 0 and 1 and a 
%  duration. The input .txt file must specify those Flash types and also 
%  one contrast value and duration for the gray interweave. 
%  For instance, if there are 4 different FlashContrasts, there are 
%  4 different 'sequences' - one Flash-to-Gray sequence for each Flash type.
% 
% Photodiode information: 
% Flash epochs are represented as light periods on the photodiode. 
% Gray interleave epochs are represented as dark periods. 
% This stimulus always starts on a flash to make the start of the first
%  epoch easy to detect from the photodiode. 
%
% Testing code for this stimulus class can be found in
%  ShortFlashOntoGrey_test3.m. Note: updates to FullFieldFlashOntoGray.m
%  are not necessarily made in ShortFlashOntoGrey_test3.m. 
% 
% Parameters:
% 	FlashDuration - durations of the non-gray interweave epochs (sec)
% 	FlashContrast - contrasts of the non-gray interweave epochs (sec)
% 	GrayDuration  - duration of gray interweave btw light or dark (sec)
%   	'flashes'
%   GrayContrast  - contrast of gray interweave 
%   RepeatRNGSeed - 0 draws a different sequence of random contrast
%   	values for the flashes, 1 draws the same random sequence each
%   	time displayStim() is called 
%
% Constructor: obj = FullFieldFlashOntoGray(txtFile)
% 	txtFile - .txt file specifying stimulus parameters and stimulus
%   	class name (see FFF_shortflash_onto_gray.txt as an example)
%
% Output (added to stim.Out):
% 	rndSeed - the integer seed used for the random number generator
%	rawStim - vector containing the stimulus contrast value played at 
%   	each frame
%   nFlashes - number of flashes presented
%   nUpdates - number of times new contrast is picked
%   stimIFI - IFI, inherited through psychtoolbox
%
% Updates:
%   3/7/17 - added definition of reconstructStim method
%   3/8/17 - added d
% last update: 3/8/17 HHY
%
%

classdef FullFieldFlashOntoGray < Stimulus 
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
        % Constructor initializes an object of this stimulus class using
        % the .txt file chosen by the user
        function obj = FullFieldFlashOntoGray(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)
            
            % Fetch parameters 
            flashDur = cell2mat(obj.FlashDuration); % vector of different flash durations
            grayDur = obj.GrayDuration{1}; % duration of gray epoch in seconds 
            flashFrames = round(flashDur./ifi);  % vector of #frames for flash types
            grayFrames = round(grayDur/ifi);       
            flashContrast = cell2mat(obj.FlashContrast); % flash contrast values 
            gray = obj.GrayContrast{1}; % grey contrast value
            nFlashType = length(flashContrast); % number of different flash types

            % Set random number generator 
            if obj.RepeatRNGSeed{1} == 1 
                % Set random number generator to produce same sequence of
                % random numbers upon each all to displayStim()
                rng('default');
            else 
                rng('shuffle'); % use a different random number sequence for each trial
            end 
            
            % Create a vector to save stimulus contrast value at each frame
            estFrames = ceil(stimDuration/ifi); % user-estimated total stimulus length
            buffer = ceil(180/ifi); % 3 min buffer
            totalFrames = estFrames + buffer; 
            rawStim = zeros(totalFrames, 1); % saves the stimulus
            nUpdates = 0; % counts number of times a new contrast is picked
            nFlashes = 0; % counts number of flashes (flash-to-grey sequences)

%             % Frame counter 
%             currFrame = 1; 
            % Counters to keep track of current position in flash-to-gray sequence 
            epochFrameCount = 0;
            isGray = 1;
            
            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window); 
            
            % MAIN LOOP over each stimulus frame
            % Keep playing the stimulus until the user presses a key on the keyboard
            for currFrame = 1:totalFrames 
                % If previous gray epoch ended, randomly draw a new flash epoch 
                if epochFrameCount == 0 && isGray == 1 
                    k = ceil(nFlashType*rand(1));
                    contrast = flashContrast(k);
                    epochFrameCount = flashFrames(k); % run new contrast value for this many frames
                    nFlashes = nFlashes + 1;
                    nUpdates = nUpdates + 1;
                    isGray = 0; % switch to flash epoch 

                % If previous flash period ended, start grey epoch
                elseif epochFrameCount == 0 && isGray == 0
                    contrast = gray;
                    epochFrameCount = grayFrames; 
                    nUpdates = nUpdates + 1;
                    isGray = 1; 
                end 

                % set contrast for photodiode  
                if (contrast == gray) % use 0.5 for gray for better visibility 
                    pdContrast = 0;
                else
                    pdContrast = 1;
                end
                    
                % Draw full field flash with the current contrast value 
                Screen('FillRect', window, contrast*obj.RGBscale, obj.screenDim); % stimulus screen
                Screen('FillRect', window, pdContrast*obj.RGBscale, obj.pdDim); % photodiode input
                % Flip the window at the framerate 
                Screen('Flip', window, startTime + currFrame*ifi);   
                pause(0)
                
                % save the contrast value at each frame
                rawStim(currFrame) = contrast;
                % update frame counter for this stimulus epoch
                epochFrameCount = epochFrameCount - 1;
                
                % If user hits keyboard, exit the stimulus loop 
                if KbCheck
                    break;
                end       
                
%                 currFrame = currFrame + 1;
            end 

            Screen('CloseAll');
            
            % save metadata specific to this class
            r = rng;
            obj.Out.rndSeed = r.Seed; % seed for random number generator 
            obj.Out.rawStim = rawStim(1:currFrame); % save stimulus contrast values to Out
            obj.Out.nFlashes = nFlashes;
            obj.Out.nUpdates = nUpdates;
            obj.Out.stimIFI = ifi;
        end 
        
        % For reconstructing the stimulus
        function [rcStim, rcStimInd] = reconstructStim(obj, ...
                stimEpochStartTimes, lightStartTimes, darkStartTimes, order)
            display('Processing photodiode data for FullFieldFlashOntoGray...')

            % lights = flashes, darks = gray
            nFlashes = length(lightStartTimes);
            nFlashTypes = length(obj.FlashContrast);
            rndSeed = obj.Out.rndSeed;
            rng(rndSeed); % seed the random number generator 
            % get sequence of FlashContrast indices
            flashInd = ceil(nFlashTypes*rand(nFlashes, 1));

            % Include gray values and sort indices based on order in time 
            % Turn the cell array of the stim class's contrast values 
            %  into an array
            flashTypes = cell2mat(obj.FlashContrast)';
            % sequence of flash contrast values
            flashArray = flashTypes(flashInd); 
            grayArray = obj.GrayContrast{1}.*...
                ones(length(darkStartTimes), 1); 

            % reconstructed stimulus as a sequence of contrast values
            contrasts = [flashArray; grayArray];
            rcStim = contrasts(order); 

            % reconstructed stimulus as a sequence of epoch indices (gray 
            %  is denoted as the highest flash index + 1
            epochInd = [flashInd; ...
                (nFlashTypes+1)*ones(length(grayArray), 1)];
            rcStimInd = epochInd(order);

            % Plot reconstructed stim over stimulus saved in the program.
            % Because saved stimulus is contrast value at each FRAME, we 
            % need to repeat each contrast value for the correct number of
            % frames depending on the duration of that flash type
    %         framesPerEpoch = floor((cell2mat(stim.obj.FlashDuration)')./stim.obj.IFI);
            framesPerEpoch = round((cell2mat(obj.FlashDuration)')./obj.IFI);
            grayFrames = round(obj.GrayDuration{1}/obj.IFI);
            % same length as flashArray
            flashRepeats = framesPerEpoch(flashInd); 
            grayRepeats = grayFrames*ones(length(darkStartTimes), 1); 
            repeatsPool = [flashRepeats; grayRepeats];
            repeats = repeatsPool(order);
            rcFrameVal = repelem(rcStim, repeats);

            % edge case - if we started with a gray period but did not
            %  detect it (would have been a black-to-dark transition on PD)
            grayContrast = obj.GrayContrast{1};
            if obj.Out.rawStim(1) == grayContrast
                rcFrameVal = repelem([grayContrast; rcStim], ...
                    [grayFrames; repeats]);
            end 

            figure;
            plot(rcFrameVal, 'LineWidth', 2.5); hold on;
            plot(obj.Out.rawStim, 'g', 'LineWidth', 1);
            legend('reconstructed stimulus', 'raw stimulus');
            ylabel('contrast'); 
            xlabel('frame number');
            title('contrast value at each stimulus frame');        
        end
        
        % selects imaging frames to be used as F0 in computing dF/F during
        %  analysis; here, last 25% of gray period
        % replaces selectBaseline_FullFieldFlashOntoGray.m
        function [baselineFrameTimes, baselineSignals] = ...
                selectDFFbaseline(obj,imgFrameTimes, bksSignal, ...
                lightStartTimes, darkStartTimes) 
            
            grayDuration = obj.GrayDuration{1};
            startBaseTime = 0.75*grayDuration;
            endBaseTime = grayDuration; 
            
            % Start times of gray interleave epochs in the stimulus.
            % Gray epochs in the stimulus are represented as 
            % dark epochs on the photodiode 
            grayStartTimes = darkStartTimes;
            
            % do not count imaging frames before first gray epoch
            % first imaging frame to count is the first imaging frame that 
            % occurs after start of gray epoch
            startImgFrame = find((imgFrameTimes - grayStartTimes(1)) > 0, ...
                1, 'first');

            % initialize arrays
            baselineFrameTimes = [];
            baselineSignals = [];

            % for each imaging frame, get time difference between 
            % transitions to gray epochs and the start of this imaging 
            %  frame 
            for i = startImgFrame:length(imgFrameTimes)
                % time difference between stimulus transitions and current 
                %  imaging frame time
                tDiffGray =  imgFrameTimes(i) - grayStartTimes;

                % only consider stimulus transitions that occured on or 
                %  before this imaging frame.
                tDiffGrayValid = tDiffGray(tDiffGray >= 0);

                % find the time from the most recent stimulus transition 
                %  to this img frame
                minGrayTimeDiff = min(tDiffGrayValid); 

                % account for empty array when full stimulus cycle hasn't 
                %  happened
                if isempty(minGrayTimeDiff)
                    minGrayTimeDiff = -1; % invalid value, evaluates to false
                end

                % if the imaging frame falls within the correct time period 
                %  (occurs after startBaseTime and before endBaseTime), 
                %  keep this imaging frame and the corresponding signal
                if ((minGrayTimeDiff > startBaseTime) && ...
                        (minGrayTimeDiff < endBaseTime))
                    % save the imaging frame time and the F value
                    baselineFrameTimes = [baselineFrameTimes; imgFrameTimes(i)];
                    baselineSignals = [baselineSignals; bksSignal(i)];
                end 

            end 
        end      
    end
end 