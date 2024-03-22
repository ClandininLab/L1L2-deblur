% SearchStimulusFlash.m
% 
% Subclass of Stimulus that presents a flashing quadrilateral in the center
%  of the screen (from the fly's perspective in degrees of visual field 
%  spanned), a set number of degrees from each edge of the screen. 
% Use to search for cells with receptive fields in the center of the
%  screen.
%
% Note: Not tested with >2 epochs. Also, assumes 1 epoch background color =
%  center color. Should probably write to be more general...
%
% Parameters:
%   Duration - duration of each epoch in seconds
%   BackgroundColor - color of edge region (i.e. background)
%   CenterColor - color of center region
%   DegFromEdge - distance from edge of center region, expressed in
%       degrees, same for all epochs
%
% Constructor: obj = SearchStimulusFlash(txtFile)
%   txtFile - .txt file specifying stimulus parameters and stimulus
%           class name
%
% Output (adds to obj.Out):
% 	rawStim - vector containing the stimulus contrast value played at 
%           each frame
%   stimIFI - stimulus IFI, inherited from psychtoolbox initialization
%
% Updates:
%   3/7/17 - added definition of reconstructStim method
%   3/8/17 - added definition of selectDFFbaseline method
% last update: 3/8/17 HHY
% 
classdef SearchStimulusFlash < Stimulus 
    properties 
        Out
        ParamList = {'Duration', 'BackgroundColor', 'CenterColor',...
            'DegFromEdge'};
        Duration
        BackgroundColor
        CenterColor
        DegFromEdge
    end 
    
    methods
        % Constructor
        function obj = SearchStimulusFlash(txtFile)
            obj@Stimulus(txtFile);
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)  
            % estimate of total number of stimulus frames based on duration
            % specified by user
            estFrames = ceil(stimDuration/ifi);
            buffer = ceil(180/ifi); % 3 min buffer
            totalFrames = estFrames + buffer;
%           % vector to save actual stimulus played.
            % column 1 is background color, column 2 is center color
            rawStim = zeros(totalFrames, 2); % preallocate 
            
            % --- compute the verticies of the quadrilateral --- %
            
            % conversion between cm and pixels
            x_cm2px = (obj.screenDim(3)-obj.screenDim(1))/...
                obj.physDist.scrnX;
            y_cm2px = (obj.screenDim(4)-obj.screenDim(2))/...
                obj.physDist.scrnY;
            
            % list of verticies
            verticies = zeros(4,2); % allocate, quadrilateral verticies

            % amount to clip from each, in pixels
            % lower left x 
            degFromEdge = obj.DegFromEdge{1}; % extract double from cell
            llx = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnX, obj.physDist.fly2scrnLL, ...
                obj.physDist.fly2scrnUL) * x_cm2px;
            % lower left y
            lly = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnY, obj.physDist.fly2scrnLL,...
                obj.physDist.fly2scrnLR) * y_cm2px;
            % lower right x
            lrx = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnX, obj.physDist.fly2scrnLR, ...
                obj.physDist.fly2scrnUR) * x_cm2px;
            % lower right y
            lry = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnY, obj.physDist.fly2scrnLR, ...
                obj.physDist.fly2scrnLL) * y_cm2px;
            % upper left x
            ulx = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnX, obj.physDist.fly2scrnUL, ...
                obj.physDist.fly2scrnUR) * x_cm2px;
            % upper left y
            uly = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnY, obj.physDist.fly2scrnUL, ...
                obj.physDist.fly2scrnUR) * y_cm2px;
            % upper right x
            urx = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnX, obj.physDist.fly2scrnUR, ...
                obj.physDist.fly2scrnLR) * x_cm2px;
            % upper right y
            ury = computeEdgeDistance(degFromEdge, ...
                obj.physDist.scrnY, obj.physDist.fly2scrnUR, ...
                obj.physDist.fly2scrnUL) * y_cm2px;
            
            % verticies, in order counterclockwise
            % lower left, lower right, upper right, upper left)
            verticies(1,:) = [obj.screenDim(1) + llx, obj.screenDim(2) + lly];
            verticies(2,:) = [obj.screenDim(1) + lrx, obj.screenDim(4) - lry];
            verticies(3,:) = [obj.screenDim(3) - urx, obj.screenDim(4) - ury];
            verticies(4,:) = [obj.screenDim(3) - ulx, obj.screenDim(2) + uly];
            
            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window);
            
            % counts total number of epochs presented (full or partial)
            epochCount = 1;
            % keep playing the stimulus until the user presses a key on the keyboard
            for frame = 1:totalFrames
                % which epoch we're currently on
                currEpoch = mod(epochCount-1,obj.nEpochs)+1;
                % update epoch counter, maybe
                epochCount = ceil((frame*ifi)/obj.Duration{currEpoch});
                
                % Draw the background onto the screen area
                Screen('FillRect',window,...
                    obj.BackgroundColor{currEpoch}*obj.RGBscale,...
                    obj.screenDim);
                % Draw the center as convex polygon defined by verticies
                Screen('FillPoly',window, ...
                    obj.CenterColor{currEpoch}*obj.RGBscale,...
                    verticies, 1);
                % Draw the input to the photodiode, does what center does
%                 Screen('FillRect', window, ...
%                     obj.CenterColor{currEpoch}*obj.RGBscale,...
%                     obj.pdDim);
                % Draw the input to the photodiode, white whenever center
                % color doesn't equal background color, black otherwise
                if obj.CenterColor{currEpoch} ~= ...
                        obj.BackgroundColor{currEpoch}
                    Screen('FillRect',window,obj.RGBscale,obj.pdDim);
                else
                    Screen('FillRect',window,0,obj.pdDim);
                end
                
                 % flip entire window at the frame rate
                Screen('Flip', window, startTime + frame*ifi);  
                
                pause(0) % super important - without it, NIDAQ won't scan
                
                % record stimulus info, background and center color
                rawStim(frame,:) = [obj.BackgroundColor{currEpoch},...
                    obj.CenterColor{currEpoch}];    
               
                % If user hits keyboard, exit the stimulus loop 
                if KbCheck
                    break;
                end    
            end 
            
            Screen('CloseAll');
            
            rawStim = rawStim(1:frame,:);
            obj.Out.rawStim = rawStim;
            obj.Out.stimIFI = ifi;
        end 
        
        % For reconstructing the stimulus
        function [rcStim, rcStimInd] = reconstructStim(obj, ...
            stimEpochStartTimes, lightStartTimes, darkStartTimes, order)
        
            display(...
                'Processing photodiode data for SearchStimulusFlash ...');

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
            contrasts = cell2mat(obj.CenterColor);
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