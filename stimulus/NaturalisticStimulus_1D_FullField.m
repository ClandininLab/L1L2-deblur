% NaturalisticStimulus_1D_FullField.m
% 
% Subclass of Stimulus that presents a full field naturalistic stimulus. 
% Note: This is written to specifically to present a different frame for
% each pattern of the 3 pattern sequence on the lightcrafter, by separating
% frames across RGB values. The encoding is specifically for the 6-bit
% setting on the lightcrafter and presents frames at 3X the lightcrafter
% frame rate.
% 
% Parameters:
% 	MatFileName - full path to .mat file containing pre-generated naturalistic
%       stimulus. File should contain one cell array with variable name
%       natStim, where each element of cell array is frame by frame
%       sequence of a single naturalistic stimulus. The number of
%       elements of the cell array is >= the number of epochs. The duration
%       of each epoch must be less than the number of frames in the
%       sequence. NOTE: I am not writing in checks into this class 
%       definition that this .mat file is correct, so ensure that the .mat
%       file is correct before running this class of stimuli
% 	Duration - duration of each epoch in seconds
%   WhichStim - which index into natStim should be presented for this epoch
%   RepeatRNGSeed - 0 draws a different sequence of random contrast
%   	values for the flashes, 1 draws the same random sequence each
%   	time displayStim() is called 
%   rgbSeq - cell array that has converted natStim into rgb frames
%
% Constructor: obj = NaturalisticStimulus_1D_FullField(txtFile) 
% 	txtFile - .txt file specifying stimulus parameters and stimulus
%   	class name
%
% Methods:
%   setParams - overrides Stimulus class's method to appropriately read in
%       Filename parameter
%   initStim - takes advantage of StimulusPerspCor abstract class's naming
%       of initStim method that is conditionally called in playStimMain,
%       though this is not a subclass of StimulusPerspCor. This is where
%       pre-loading of the stimulus happens, before anything is displayed.
%       Converts intensity sequence in natStim to RGB frames. Load into 
%
% Output (adds to stim.Out):
% 	rawStim - vector containing rgb frame contrast value played at 
%       each Psychtoolbox frame
% 	rndSeed - the integer seed used for the random number generator
%   numEpochsPresented - number of natural stimulus epochs presented
%   epochsPresented - sequence of which epochs presented
%   stimIFI - IFI of stimulus presentation, inherited from psychtoolbox
%   Missed - psychtoolbox readout of missed frames (for debugging)
%   VBLTimestamp - psychtoolbox VBL time stamps for each frame (for
%       debugging)
%   StimulusOnsetTime - psychtoolbox output of StimulusOnsetTime (for
%       debugging)
%
% Updates:
% 7/12/19 HHY
% last update: 8/9/19 HHY - fixed some minor bugs
% 
classdef NaturalisticStimulus_1D_FullField < Stimulus 
    properties 
        Out
        ParamList = {'Duration', 'WhichStim', 'RepeatRNGSeed'}; 
        Duration
        WhichStim
        RepeatRNGSeed
        MatFileName
        rgbSeq
    end 
    methods
        % Constructor initializes an object of this stimulus class using
        % the .txt file chosen by the user
        function obj = NaturalisticStimulus_1D_FullField(txtFile)
            obj@Stimulus(txtFile);
        end
        
        % reads in .txt file and sets parameters. Overwrites one from
        %  Stimulus class
        function setParams(obj)
            fileID = fopen(obj.TxtFileName); % open the txt file
            line1 = textscan(fileID, '%s %f', 1); % number of params
            line2 = textscan(fileID, '%s %f', 1); % number of epochs
            line3 = textscan(fileID, '%s %s', 1); % stimulus function name
            line4 = textscan(fileID, '%s %s', 1); % filename
           
            % save the values
            obj.nParams = line1{2};  
            obj.nEpochs = line2{2};   
%             className = line3{2};
            
%             % check to ensure stimulus class in text file matches this one
%             if ~strcmp(className, 'FullFieldFlash')
%                 fprintf('Error! Class name mismatch');
%                 return;
%             end 
            obj.MatFileName = line4{2}{1};
            
            % loop through txt file and list of parameters to assign values 
            % to the expected parameters 
            while ~feof(fileID)
                pstr = textscan(fileID,'%s', 1); % name of parameter
                paramName = char(pstr{1}); % access the cell to save the name 
                formatStr = repmat('%f ', 1, obj.nEpochs); % assumes user set "EPOCHS" correctly
                values = textscan(fileID, formatStr, 1); % parameter values
                
                % loop through parameter list
                for i = 1:length(obj.ParamList) 
                    % if parameter name in txt file matches name in class's
                    % parameter list, assign the values of the paramters to
                    % the corresponding property in this class
                    if strcmp(paramName, obj.ParamList{i})
                        % get rid of empty cell elements and NaNs
                        pcell = values(1:end);
                        pcell = pcell(~cellfun('isempty', pcell));
                        pmat = cell2mat(pcell);
                        pcell = pcell(~isnan(pmat)); 
                        % save the non-Nan param values
                        pname = matlab.lang.makeValidName(paramName); 
                        obj.(pname) = pcell;
                    end 
                end 
            end 
            fclose(fileID);
        end
        
        % initializes stimulus by loading sequence of intensity values from
        %   file, converts them to RGB frames
        function initStim(obj,window,ifi)
            load(obj.MatFileName, 'natStim');
            
            for i = 1:obj.nEpochs
                obj.rgbSeq{i} = convertInt2RGB(natStim{obj.WhichStim{i}});
            end
        end
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)  
            
            % Seconds per stimulus epochs, as array
            epochDurs = cell2mat(obj.Duration);
            % number of frames per epoch
            epochFrames = round(epochDurs./ifi);
            % number of epochs
            numEpochs = length(epochDurs);
            
            % Set random number generator 
            if obj.RepeatRNGSeed{1} == 1 
                % Set random number generator to produce same sequence of
                % random numbers upon each all to displayStim()
                rng('default');
            else 
                % use a different random number sequence for each trial
                rng('shuffle'); 
            end 
            
            % estimate of total number of stimulus frames based on duration
            % specified by user
            estFrames = ceil(stimDuration/ifi);
            buffer = ceil(180/ifi); % 3min buffer
            totalFrames = estFrames + buffer;
            % vector to save actual stimulus played
            rawStim = zeros(totalFrames, 3); 
            
            % Counters to keep track of current frame in epoch, which epoch 
            epochFrameCount = 1;
            % pick random starting epoch
            whichEpoch = ceil(numEpochs*rand(1)); 
            % number of epoch presented (doesn't have to be complete)
            numEpochPresented = 1;
            % list of epochs presented
            epochsPresented = whichEpoch;
            
            % Arrays for saving outputs of Screen('Flip')
            missed = zeros(totalFrames, 1);
            vblTimes = zeros(totalFrames, 1);
            stimOnsetTimes = zeros(totalFrames, 1);

            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window);
            
            % keep playing the stimulus until the user presses a key on the keyboard
            for currFrame = 1:totalFrames 
                
                % if previous epoch is finished
                if (epochFrameCount > epochFrames(whichEpoch))
                    % pick a new epoch
                    whichEpoch = ceil(numEpochs*rand(1));
                    % update counter of number of epochs presented
                    numEpochPresented = numEpochPresented + 1;
                    % update array keeping track of which epochs presented
                    epochsPresented = cat(2, epochsPresented, whichEpoch);
                    % reset epoch frame counter
                    epochFrameCount = 1;
                end

                % get current stimulus frame
                currStimFrame = ...
                    obj.rgbSeq{whichEpoch}(epochFrameCount,:);
                    
                % change photodiode signal depending on which frame of
                %  epoch
                switch epochFrameCount
                    case 1 % 1st frame
                        % photodiode signal for 1st frame is 111
                        currPDFrame = convertInt2RGB([1 1 1]);
                    case 2 % 2nd frame
                        % photodiode signal for 2nd frame is 101
                        currPDFrame = convertInt2RGB([1 0 1]);
                    otherwise % all other frames
                        % photodiode signal for all other frames is 100
                        currPDFrame = convertInt2RGB([1 0 0]);
                end

                % advance epoch frame counter
                epochFrameCount = epochFrameCount + 1;
                
                % Draw the stimulus 
                Screen('FillRect', window, currStimFrame, obj.screenDim); 
                % Draw the input to the photodiode 
                Screen('FillRect', window, currPDFrame, obj.pdDim);
                
                 % flip entire window at the frame rate
                [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed]...
                    = Screen('Flip', window, startTime + currFrame*ifi);  
                
                pause(0) % super important - without it, NIDAQ won't scan
                
                % save current frame
                rawStim(currFrame,:) = currStimFrame; 
                
                % save PTB timing information
                missed(currFrame) = Missed;
                vblTimes(currFrame) = VBLTimestamp;
                stimOnsetTimes(currFrame) = StimulusOnsetTime;
                
                % If user hits keyboard, exit the stimulus loop 
                if KbCheck
                    break;
                end    
            end 
            
            Screen('CloseAll');
            
            
            % save stimulus-specific output variables
            rawStim = rawStim(1:currFrame,:);
            r = rng;
            
            obj.Out.rawStim = rawStim;
            obj.Out.rndSeed = r.Seed;
            obj.Out.numEpochPresented = numEpochPresented;
            obj.Out.epochsPresented = epochsPresented;
            
            obj.Out.stimIFI = ifi;
            
            obj.Out.Missed = missed(1:currFrame);
            obj.Out.VBLTimestamp = vblTimes(1:currFrame);
            obj.Out.StimulusOnsetTime = stimOnsetTimes(1:currFrame);
        end 
        
        % For reconstructing the stimulus
        % rcStimInd returns start times of epoch, instead of indicies
        % rcStim as cell array of intensity values for each epoch
        function [rcStim, rcStimInd] = reconstructStim(obj, ...
        	stimEpochStartTimes, lightStartTimes, darkStartTimes, order)
            
            display('Processing photodiode data for NaturalisticStimulus_1D_FullField ...');

            % epoch start is when precisely 4 light patterns are presented
            % sequentially
            % initialize array for epoch starts
            epochStartTimes = [];
            
            % check if each light start is an epoch start
            for i = 1:length(lightStartTimes)
                % 3 light starts within 3 patterns, with 1/2 pattern
                %  buffer
                chkStart = lightStartTimes(i);
                chkEnd = lightStartTimes(i) + (obj.Out.stimIFI * (7/6));
                numLightStarts3 = length(find((lightStartTimes > chkStart) ...
                    & (lightStartTimes <= chkEnd)));
                
                % 3 light starts within 4 patterns, with 1/2 pattern buffer
                chkEnd = lightStartTimes(i) + (obj.Out.stimIFI * (9/6));
                numLightStarts4 = length(find((lightStartTimes > chkStart) ...
                    & (lightStartTimes <= chkEnd)));
                
                if ((numLightStarts3 == 3)&& (numLightStarts4 == 3))
                    epochStartTimes = [epochStartTimes; lightStartTimes(i)];
                end
            end
            
            % assign rcStimInd to be epoch start times just for this
            %  stimulus
            rcStimInd = epochStartTimes;
            
            % reconstruct stimulus presented for each epoch, except last
            %  epoch, which is partial
            
            % duration of each epoch, in frames
            epochDurs = round(diff(epochStartTimes) / obj.Out.stimIFI);
            epochDurs(epochDurs>1000) = 1000; %220518 MMP
            
            % preallocate cell array for reconstructed sequence
            recSeq = cell(length(epochStartTimes),1);
            
            for i = 1:length(epochDurs)
                % which epoch type
                whichEpochType = obj.Out.epochsPresented(i);
                
                % RGB sequence for this epoch, of appropriate duration
                rgbSeqEpoch = obj.rgbSeq{whichEpochType}(1:epochDurs(i),:);
                
                % convert to intensity sequence & add to cell array
                recSeq{i} = convertRGB2Int(rgbSeqEpoch);  
            end
            
            % reconstruct last epoch
            % end time of last epoch as last light start + 1 frame
            lastLightStartTime = lightStartTimes(end);
            
            % duration of last epoch, in frames 
            lastEpochDur = round(lastLightStartTime - ...
                epochStartTimes(end)) + 1;
            
            % RGB sequence for this epoch
            rgbSeqEpoch = obj.rgbSeq{obj.Out.epochsPresented(...
                length(epochDurs) + 1)}(1:lastEpochDur,:);
            
            % convert to intensity sequence and add to cell array
            recSeq{end} = convertRGB2Int(rgbSeqEpoch);
            
            % assign rcStim as reconstructed stimulus cell array
            rcStim = recSeq;
            
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
