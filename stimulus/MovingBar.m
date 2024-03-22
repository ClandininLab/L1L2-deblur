% MovingBar.m
% 
% Stimulus that presents a single bar that moves orthogonal to its long
%   axis.
% Perspective corrected to account for angular distortion of screen
%   position. Angular extent in all directions is corrected for. See
%   PerspectiveCorrectionTest_v1.m for example of implementation of
%   perspective correction.
% Each epoch is single bout of the bar moving, with defined parameters.
% Order of epochs is randomized without replacement. That order repeats
%  until stimulus is stopped.
% 
% Photodiode information: 
%   Each epoch begins with 2 flashes. Otherwise, black.
% Note: presentation of first epoch start is always off, as visible from
%   photodiode signal; don't use first bout in data analysis
%   update: maybe not after restructuring with initStim
% 
% Parameters:
%   BackgroundContrast - contrast of background (0-1), for each epoch
%   BarContrast - contrast of bar (0-1), for each epoch
%   BarWidth - width of bar (in degrees of fly's visual field), for each
%       epoch
%   BarAngle - angle of bar (in degrees where 0 deg is vertical, and
%       positive angles rotate CCW); not meaningful to exceed 180 deg,
%       which is equivalent to 0 deg
%   BarVelocity - velocity of bar (in degrees of fly's visual field per
%       second); negative velocities move bar in opposite direction of
%       positive
%   EpochDuration - how long (in seconds) bout is presented for
%
% Input: (see movingBar_darkC1_5deg_4dir_20ds.txt as an example)
%   txtFile - .txt file specifying stimulus parameters and stimulus
%   	class name
%
% Output (stim.Out):
%   degShift - number of degrees texture shifted, played at each frame
%   whichEpoch - which epoch was presented at each frame
%
% Updates:
%   1/15/17 HHY
%   3/28/17 - HHY - add resconstructStim and selectDFFbaseline functions
%

classdef MovingBar < StimulusPerspCor 
    properties 
        Out
        Init
        ParamList = {'BackgroundContrast', 'BarContrast', 'BarWidth',...
            'BarAngle', 'BarVelocity', 'EpochDuration'};
        BackgroundContrast 
        BarContrast
        BarWidth
        BarAngle
        BarVelocity
        EpochDuration

    end 
    methods
        % Constructor initializes an object of this stimulus class using
        % the .txt file chosen by the user
        function obj = MovingBar(txtFile)
            obj@StimulusPerspCor(txtFile);
        end 
        
        % Stimulus initialization code, pre-generates needed info for
        % displaying stimulus
        function initStim(obj,window,ifi)
            % -- Fetch parameters - convert to vector, 1 value per epoch
            bkgrdContrast = cell2mat(obj.BackgroundContrast);
            barContrast = cell2mat(obj.BarContrast);
            barWidth = cell2mat(obj.BarWidth);
            barAngle = cell2mat(obj.BarAngle);
            barVelocity = cell2mat(obj.BarVelocity);
            epochDur = cell2mat(obj.EpochDuration);
            % number of frames for each epoch
            epochFrames = round(epochDur./ifi); 
            
            % -- Compute static measurements for mapping texture to screen
            
            % compute vertices of persepctive quadrilateral
            [urx, ury, ulx, uly, llx, lly, lrx, lry] = ...
                computePerspecQuadVerticies(obj.physDist.fly2scrnUL, ...
                obj.physDist.fly2scrnUR, obj.physDist.fly2scrnLL, ...
                obj.physDist.fly2scrnLR, obj.physDist.scrnX, ...
                obj.physDist.scrnY, obj.physDist.fly2horizL, ...
                obj.physDist.fly2horizR, obj.physDist.fly2topL,...
                obj.physDist.fly2topR, obj.txtrParams.xCorner,...
                obj.txtrParams.yCorner);
            
            % convert verticies of quadrilateral to coefficients for 
            %  mapping quadrilateral to unit square
            xVert = [llx; lrx; urx; ulx];
            yVert = [lly; lry; ury; uly];
            [xCoeffs, yCoeffs] = mapPerspecQuadVertToUnitSquare(...
                xVert,yVert);
            
            % Use the coefficients to get (x,y) coordinates for getting  
            %  intensity values of texture
            screenXDim = obj.screenDim(3) - obj.screenDim(1);
            screenYDim = obj.screenDim(4) - obj.screenDim(2);

            % get (l,m) at resolution of txtrRes
            numBinsX = screenXDim / obj.txtrParams.txtrRes;
            numBinsY = screenYDim / obj.txtrParams.txtrRes;
            
            % compute l and m: evenly spaced between 0 and 1
            [m,l] = meshgrid((0:(1/(numBinsY-1)):1),...
                (0:(1/(numBinsX-1)):1));
    
            % get (x,y) corresponding to (l,m)
            [quadX,quadY] = LMtoXY(xCoeffs,yCoeffs,l,m);
            
            % coordinates of test texture 
            [txtrX, txtrY] = meshgrid(...
                (1/obj.txtrParams.txtrRes):(1/obj.txtrParams.txtrRes):...
                obj.txtrParams.txtrSize,...
                (1/obj.txtrParams.txtrRes):(1/obj.txtrParams.txtrRes):...
                obj.txtrParams.txtrSize);
    
            
            % -- Pre-generate all textures to be presented
            % will hold pointers to textures for all epochs
            txtrPointers = cell(obj.nEpochs,1);
            % will hold all shifts, by epoch
            epochDegShift = cell(obj.nEpochs,1);
            for i = 1:obj.nEpochs
                
                % -- create starting texture --
                
                % minimum starting texture size, in degrees, to account for 
                %  rotation
                minRotTxtrSize = ceil(2*obj.txtrParams.txtrSize*sind(45));

                % bar width must be whole number of pixels
                barWidthPx = round(barWidth(i)*obj.txtrParams.txtrRes);
                
                % minimum starting texture size, to account for duration of
                %  epoch
                minTimeTxtrSize = barWidthPx/obj.txtrParams.txtrRes + ...
                    abs(barVelocity(i))*epochDur(i);                
                
                % starting texture size is larger of 2 minimums; convert to
                %  pixels from degrees
                startTxtrSize = max([minRotTxtrSize, minTimeTxtrSize])*...
                    obj.txtrParams.txtrRes;
                % starting texture; define bar in row dimension; starts
                % from (0,0) and extends barWidth * txtrRes rows by
                % startTxtrSize columns
                startTxtr = [ones(barWidthPx,startTxtrSize)...
                    *barContrast(i); ...
                    ones(startTxtrSize-barWidthPx,...
                    startTxtrSize)*bkgrdContrast(i)];
                
                % -- generate sequence of textures --
                
                % preallocate for texture pointers for this epoch
                thisTxtrPoint = zeros(epochFrames(i),1);
                % preallocate for degShift for this epoch
                thisDegShift = zeros(epochFrames(i),1);
                for j = 1:epochFrames(i)
                    % amount to shift on each frame, round to nearest pixel
                    shiftAmount = round((j-1)*barVelocity(i)*ifi*...
                        obj.txtrParams.txtrRes);
                    thisDegShift(j) = shiftAmount/obj.txtrParams.txtrRes;
                    % shift texture
                    shiftTxtr = circshift(startTxtr,[shiftAmount 0]);
                    % rotate texture, only if there is rotation to be done
                    if barAngle(i)
                        rotTxtr = imrotate(shiftTxtr,barAngle(i),...
                            'bilinear');
                    else
                        rotTxtr = shiftTxtr;
                    end
                    % crop texture to txtrSize
                    rotTxtrSize = size(rotTxtr);
                    startX = floor(rotTxtrSize(1)/2 - ...
                        (obj.txtrParams.txtrSize*obj.txtrParams.txtrRes)/2);
                    endX = startX + ...
                        (obj.txtrParams.txtrSize*obj.txtrParams.txtrRes) - 1;
                    startY = floor(rotTxtrSize(2)/2 - ...
                        (obj.txtrParams.txtrSize*obj.txtrParams.txtrRes)/2);
                    endY = startY + ...
                        (obj.txtrParams.txtrSize*obj.txtrParams.txtrRes) - 1;
                    txtr = rotTxtr(startX:endX,startY:endY);
                    
                    % convert coordiates and texture to perspective 
                    %  corrected texture
                    pcTxtr = interp2(txtrX,txtrY,txtr,quadX,quadY);
                    
                    % convert perspective corrected texture into PTB
                    %  texture
                    thisTxtrPoint(j) = Screen('MakeTexture',window,...
                        pcTxtr*obj.RGBscale);
                end
                % save all texture pointers
                txtrPointers{i} = thisTxtrPoint;
                % save all shifts
                epochDegShift{i} = thisDegShift;
            end
            % save to object
            obj.Init.txtrPointers = txtrPointers;
            obj.Init.epochDegShift = epochDegShift;
        end
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi, stimDuration)
            % -- Fetch parameters - convert to vector, 1 value per epoch
            epochDur = cell2mat(obj.EpochDuration);
            % number of frames for each epoch
            epochFrames = round(epochDur./ifi); 
            % randomize order of epochs
            epochOrder = randperm(obj.nEpochs);
            
            
            % -- Create a vector to save --
            % time to run through all epochs once
            totEpochTime = sum(epochDur);
            % number of times to run through all epochs, based on
            %  user-estimated total stimulus length, rounding up
            numReps = ceil(stimDuration / totEpochTime);
            % total number of frames to present 
            totalFrames = ceil(numReps*totEpochTime/ifi);
            
            % pre-allocate to save degShift and which epoch
            degShift = zeros(totalFrames,1);
            whichEpoch = zeros(totalFrames,1);

            % Frame counter 
            frameCount = 1; 
            % flag to break out of presentation
            endFlag = 0;
            
            % --- Psychtoolbox code --- %
            startTime = Screen('Flip', window);
            
            % loop over number of reps of whole stimulus
            for i = 1:numReps
                % loop over each epoch
                for j = 1:obj.nEpochs
                    curEpoch = epochOrder(j);
                    % loop through frames of each epoch
                    for k = 1:epochFrames(curEpoch)
                        % define photodiode contrast - white for first 2
                        %  frames of epoch, black otherwise
                        if (k<3)
                            pdContrast = 1;
                        else
                            pdContrast = 0;
                        end
                        
                        % draw stimulus
                        % texture
                        Screen('DrawTexture',window,...
                            obj.Init.txtrPointers{curEpoch}(k),[],...
                            obj.screenDim,[],0);
                        % photodiode 
                        Screen('FillRect', window, ...
                            pdContrast*obj.RGBscale, obj.pdDim); 
                        
                        % Flip the window at the framerate 
                        Screen('Flip', window, startTime + frameCount*ifi);   
                        pause(0)
                        
                        % save current values of degShift and epoch
                        degShift(frameCount) = ...
                            obj.Init.epochDegShift{curEpoch}(k);
                        whichEpoch(frameCount) = curEpoch;
                        
                        % If user hits keyboard, exit the stimulus loop 
                        if KbCheck
                            endFlag = 1;
                            break;
                        end
                        
                        % update frame counter
                        frameCount = frameCount + 1;
                    end
                    % user hit keyboard, exit stimulus loop
                    if endFlag
                        break;
                    end
                end
                % user hit keyboard, exit stimulus loop
                if endFlag
                    break;
                end
            end
            
            Screen('CloseAll')
            
            % save metadata
            if endFlag % end early, save only used portion
                obj.Out.degShift = degShift(1:(frameCount-1));
                obj.Out.whichEpoch = whichEpoch(1:(frameCount-1));
            else
                obj.Out.degShift = degShift;
                obj.Out.whichEpoch = whichEpoch;          
            end
        end 
        
        % Reconstructs the stimulus presented, for analysis
        %  adds epochTransTimes and epochsPresented to obj.Out
        %  rcStim and rcStimInd are the same and both reference the epoch
        function [rcStim, rcStimInd] = reconstructStim(obj, ...
                stimEpochStartTimes, lightStartTimes, darkStartTimes, ...
                order)
            display('Processing photodiode data for MovingBar...')
            
%             % epoch changes only at light preceeded by long dark
%             diffLightStartTimes = diff(lightStartTimes);
%             
%             % ignore second flash presented at start of epoch
%             % get indicies of lightStartTimes where epoch really changed
%             transInd = find(diffLightStartTimes > (2*obj.IFI))+1;
%             transInd = [1; transInd];
%             
%             % epoch transition times
%             epochTransTimes = lightStartTimes(transInd);

            % epoch transitions at start of light flash
            epochTransTimes = lightStartTimes;
            
            % get epochs presented (which epoch is starting)
            if (obj.Out.whichEpoch == 0)
                epochInd = [1; (find(diff(obj.Out.whichEpoch(1:(end-1)))) + 1)];
                epochs = obj.Out.whichEpoch(epochInd);
            else
                epochInd = [1; (find(diff(obj.Out.whichEpoch)) + 1)];
                epochs = obj.Out.whichEpoch(epochInd);
            end
                
            % save epoch transition times, epochs presented
            obj.Out.epochTransTimes = epochTransTimes;
            obj.Out.epochsPresented = epochs;
            
            % compute rcStim and rcStimInd
            %  not actually that useful for this stimulus, but used for
            %  prelim plotting
            rcStimInd = zeros(size(stimEpochStartTimes));
            rcStim = zeros(size(stimEpochStartTimes));
            for i=1:length(stimEpochStartTimes)
                currTime = stimEpochStartTimes(i);
                
                % find nearest epoch transition that happened at the same
                %  time or before currTime
                whichInd = find((currTime >= epochTransTimes),1,'last');
                % use to index into epochs to get which epoch it is on
                rcStimInd(i) = epochs(whichInd);
                % normalize to max of 1
                rcStim(i) = rcStimInd(i)/max(epochs);
            end
            
        end
        
        % selects imaging frames to be used as F0 in computing dF/F during
        %  analysis
        % for this stimulus, first baselineFraction (25%) of each epoch
        function [baselineFrameTimes, baselineSignals] = ...
                selectDFFbaseline(obj,imgFrameTimes, bksSignal, ...
                lightStartTimes, darkStartTimes)
            % what fraction of epoch to use for baseline (from start) 
            baselineFraction = 0.25;
            
            % initialize output variables
            baselineFrameTimes = [];
            baselineSignals = [];
            
            % get epoch durations as array
            epochDur = cell2mat(obj.EpochDuration);
            
            % start counting imaging frames starting for those that fall
            %  after start of first epoch
            startImgFrame = find((imgFrameTimes - ...
                obj.Out.epochTransTimes(1)) > 0, 1, 'first');
            
            for i=startImgFrame:length(imgFrameTimes)
                currTime = imgFrameTimes(i);
                
                % find nearest epoch transition that happened at the same
                %  time or before currTime
                whichInd = find((currTime >= obj.Out.epochTransTimes),...
                    1,'last');
                % use to index into obj.Out.epochsPresented to get which
                %  epoch that imaging frame fell into
                whichEpoch = obj.Out.epochsPresented(whichInd);
                
                % define start and end time points for baseline frames
                currEpochDur = epochDur(whichEpoch);
                endBaseTime = currEpochDur * baselineFraction;
                startBaseTime = 0; % start always at 0, first baselineFraction
                
                % time relative to epoch start
                relTime = currTime - obj.Out.epochTransTimes(whichInd);
                
                % is this frame in first baselineFraction of this epoch?
                %  i.e. occurs after startTime and before endTime
                if ((relTime >= startBaseTime) && (relTime < endBaseTime))
                    % save the imaging frame time and the F value
                    baselineFrameTimes = [baselineFrameTimes; ...
                        imgFrameTimes(i)];
                    baselineSignals = [baselineSignals; bksSignal(i)];
                end
            end
        end
    end 
end 