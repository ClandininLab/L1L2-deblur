% StaticAndMovingSquareGrating.m
% 
% Stimulus that presents gray, a static square grating, and then the same
%   grating moving orthogonal to the bars.
% Perspective corrected to account for angular distortion of screen
%   position. Angular extent in all directions is corrected for. See
%   PerspectiveCorrectionTest_v1.m for example of implementation of
%   perspective correction.
% Each epoch is gray, static, moving, with defined parameters.
% Order of epochs is randomized without replacement. That order repeats
%  until stimulus is stopped.
% 
% Photodiode information: 
%   Gray begins with 2 flashes, static with 1, moving with 3. Otherwise,
%   black.
% 
% Parameters:
%   DarkContrast - contrast of dark bars (0-1), for each epoch
%   LightContrast - contrast of light bars (0-1), for each epoch
%   GrayContrast - contrast of gray period (0-1), for each epoch
%   GratingPeriod - period of square grating, in degrees; equal to 2X the
%       width of a bar
%   GratingAngle - angle of grating (in degrees where 0 is vertical, and
%       positive angles rotate CCW). Avoid exceeding 180 deg, as this will
%       offset the grating by 1/2 period but not change the direction of
%       movement.
%   GratingVelocity - velocity of grating during movement period (in
%       degrees of fly's visual field per second); negative velocities move
%       bar in opposite direction of positive
%   PhaseOffset - initial offset in position of grating, in degrees of
%       fly's visual field
%   GrayDuration - duration (in seconds) of the gray period, for each epoch
%   StaticDuration - duration (in seconds) of the static grating
%       presentation period, for each epoch
%   MovingDuration - duration (in seconds) of the moving grating
%       presentaiton period, for each epoch
%
%
% Input: (see staticMovingGrating_20deg_20ds_c1_12dir.txt as an example)
%   txtFile - .txt file specifying stimulus parameters and stimulus
%   	class name
%
% Output (stim.Out):
%   degShift - number of degrees texture shifted, played at each frame.
%       Will be 0 when stimulus presented is gray or static.
%   whichStim - which stimulus (gray, static, moving = 1, 2, 3) was
%       presented at each frame
%   whichEpoch - which epoch was presented at each frame
%
% Updates:
%   1/16/17 - HHY
%   3/28/17 - HHY - add resconstructStim and selectDFFbaseline functions
%

classdef StaticAndMovingSquareGrating < StimulusPerspCor 
    properties 
        Out
        Init
        ParamList = {'DarkContrast', 'LightContrast', 'GrayContrast',...
            'GratingPeriod', 'GratingAngle', 'GratingVelocity',...
            'PhaseOffset','GrayDuration','StaticDuration',...
            'MovingDuration'};
        DarkContrast
        LightContrast
        GrayContrast
        GratingPeriod
        GratingAngle
        GratingVelocity
        PhaseOffset
        GrayDuration
        StaticDuration
        MovingDuration

    end 
    methods
        % Constructor initializes an object of this stimulus class using
        % the .txt file chosen by the user
        function obj = StaticAndMovingSquareGrating(txtFile)
            obj@StimulusPerspCor(txtFile);
        end 
        
        % Stimulus initialization code, pre-generates needed info for
        % displaying stimulus
        function initStim(obj,window,ifi)
            % -- Fetch parameters - convert to vector, 1 value per epoch
            darkContrast = cell2mat(obj.DarkContrast);
            lightContrast = cell2mat(obj.LightContrast);
            gratingPeriod = cell2mat(obj.GratingPeriod);
            gratingAngle = cell2mat(obj.GratingAngle);
            gratingVelocity = cell2mat(obj.GratingVelocity);
            phaseOffset = cell2mat(obj.PhaseOffset);
            movingDuration = cell2mat(obj.MovingDuration);
            % number of frames of moving grating
            movingNumFrames = round(movingDuration./ifi);

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
            % pointers for textures
            txtrPointers = cell(obj.nEpochs,1);
            % will hold all shifts, by epoch; for moving only
            epochDegShift = cell(obj.nEpochs,1);
            for i = 1:obj.nEpochs
                
                % -- create starting texture --
                
                % minimum starting texture size, in degrees, to account for 
                %  rotation
                minRotTxtrSize = ceil(2*obj.txtrParams.txtrSize*sind(45));
                % to wrap grating correctly, texture must span integer
                %  number of periods
                nPeriods = ceil(minRotTxtrSize / gratingPeriod(i));
                % to generate texture, gratingPeriod/2 * texture resolution
                %  must be integer. Round to nearest integer. Width of
                %  single bar in pixels
                barWidth = round(gratingPeriod(i)/2 * ...
                    obj.txtrParams.txtrRes);
                % convert to starting texture size, in pixels
                startTxtrSize = nPeriods * barWidth*2;
                
                % starting texture; define grating in row dimension; starts
                % from (0,0) with black then white as 1 period
                txtrOnePeriod = ...
                    [ones(barWidth,startTxtrSize) * darkContrast(i); ...
                    ones(barWidth,startTxtrSize) * lightContrast(i)];
                startTxtr = repmat(txtrOnePeriod,nPeriods,1);
                
                % if there's a phase offset, shift the starting texture
                %  appropriately
                if phaseOffset(i)
                    phaseShiftAmount = round(phaseOffset(i)*...
                        obj.txtrParams.txtrRes);
                    startTxtr = circshift(startTxtr,[phaseShiftAmount 0]);
                end
                
                % -- generate sequence of textures --
                % amount to shift on each frame, in pixels
                shiftAmount = 0;
                % period, in pixels
                periodPx = barWidth * 2;
                frameCount = 1;
                totalShifted = 0;

                % preallocate for texture pointers for this epoch
                thisTxtrPoint = zeros(movingNumFrames(i),1);
                % preallocate for degShift for this epoch
                thisDegShift = zeros(movingNumFrames(i),1);
                while (mod(shiftAmount,periodPx)||(totalShifted==0))
                    % amount to shift on each frame, round to nearest pixel
                    shiftAmount = round((frameCount-1)*gratingVelocity(i)*...
                        ifi*obj.txtrParams.txtrRes);
                    % amount shifted, in degrees
                    thisDegShift(frameCount) = shiftAmount/...
                        obj.txtrParams.txtrRes;
                    % shift texture
                    shiftTxtr = circshift(startTxtr,[shiftAmount 0]);
                    % rotate texture, only if there is rotation to be done
                    if gratingAngle(i)
                        rotTxtr = imrotate(shiftTxtr,gratingAngle(i),...
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
                    thisTxtrPoint(frameCount) = Screen('MakeTexture',...
                        window,pcTxtr*obj.RGBscale);
                    
                    % update frame counter
                    frameCount = frameCount + 1;
                    % update shift counter
                    totalShifted = shiftAmount;
                end
                % crop texture pointers and epoch shifts vector; subtract 2
                %  1 for last update of frameCount, 1 for fact that first
                %  and last textures are the same
                thisTxtrPoint = thisTxtrPoint(1:(frameCount-2));
                thisDegShift = thisDegShift(1:(frameCount-2));
                
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
            grayContrast = cell2mat(obj.GrayContrast);
            grayDuration = cell2mat(obj.GrayDuration);
            % number of frames of gray
            grayNumFrames = round(grayDuration./ifi);
            staticDuration = cell2mat(obj.StaticDuration);
            % number of frames of static grating
            staticNumFrames = round(staticDuration./ifi);
            movingDuration = cell2mat(obj.MovingDuration);
            % number of frames of moving grating
            movingNumFrames = round(movingDuration./ifi);
            
            % randomize order of epochs
            epochOrder = randperm(obj.nEpochs);
            
            
            % -- Create a vector to save --
            % time to run through all epochs once
            totEpochTime = sum(grayDuration)+sum(staticDuration)+...
                sum(movingDuration);
            % number of times to run through all epochs, based on
            %  user-estimated total stimulus length, rounding up
            numReps = ceil(stimDuration / totEpochTime);
            % total number of frames to present 
            totalFrames = ceil(numReps*totEpochTime/ifi);
            
            % pre-allocate to save degShift and which epoch
            degShift = zeros(totalFrames,1);
            whichStim = zeros(totalFrames,1);
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
                    
                    % gray period of epoch
                    for k = 1:grayNumFrames(curEpoch)
                        % If user hits keyboard, exit the stimulus loop; 
                        %  also end when other loops broken out of
                        if (KbCheck || endFlag)
                            endFlag = 1;
                            break;
                        end
                        
                        % define photodiode contrast - white for first 2
                        %  frames of epoch, black otherwise
                        if (k<3)
                            pdContrast = 1;
                        else
                            pdContrast = 0;
                        end
                        
                        % draw gray screen
                        Screen('FillRect', window, ...
                            grayContrast(curEpoch)*obj.RGBscale, ...
                            obj.screenDim);
                        % photodiode 
                        Screen('FillRect', window, ...
                            pdContrast*obj.RGBscale, obj.pdDim); 
                        
                        % Flip the window at the framerate 
                        Screen('Flip', window, startTime + frameCount*ifi);   
                        pause(0)
                        
                        % save current values of degShift and epoch
                        degShift(frameCount) = 0; % gray, no shift
                        whichStim(frameCount) = 1; % gray is 1
                        whichEpoch(frameCount) = curEpoch;
                        
                        % update frame counter
                        frameCount = frameCount + 1;
                    end
                    
                    % static grating period of epoch
                    for k = 1:staticNumFrames(curEpoch)
                        % If user hits keyboard, exit the stimulus loop; 
                        %  also end when other loops broken out of
                        if (KbCheck || endFlag)
                            endFlag = 1;
                            break;
                        end
                        
                        % define photodiode contrast - white for first
                        %  frame of epoch, black otherwise
                        if (k<2)
                            pdContrast = 1;
                        else
                            pdContrast = 0;
                        end
                        
                        % draw static texture; always first texture for
                        %  that epoch
                        Screen('DrawTexture',window,...
                            obj.Init.txtrPointers{curEpoch}(1),[],...
                            obj.screenDim,[],0);
                        
                        % photodiode 
                        Screen('FillRect', window, ...
                            pdContrast*obj.RGBscale, obj.pdDim);   
                        
                        % Flip the window at the framerate 
                        Screen('Flip', window, startTime + frameCount*ifi);   
                        pause(0)
                        
                        % save current values of degShift and epoch
                        degShift(frameCount) = 0; % static, no shift
                        whichStim(frameCount) = 2; % static is 2
                        whichEpoch(frameCount) = curEpoch;
                        
                        % update frame counter
                        frameCount = frameCount + 1;                                                
                    end
                    
                    % moving grating period of epoch
                    for k = 1:movingNumFrames(curEpoch)
                        % If user hits keyboard, exit the stimulus loop; 
                        %  also end when other loops broken out of
                        if (KbCheck || endFlag)
                            endFlag = 1;
                            break;
                        end
                        
                        % define photodiode contrast - white for first 3
                        %  frames of epoch, black otherwise
                        if (k<4)
                            pdContrast = 1;
                        else
                            pdContrast = 0;
                        end
                        
                        % which texture - 1st frame presented is 2nd 
                        %  texture, which follows correctly from static;
                        %  wraps around dealing w/indexing appropriately
                        curTxtrInd = mod(k,...
                            length(obj.Init.txtrPointers{curEpoch}))+1;
                        % draw moving texture; loop through available
                        %  textures for that epoch
                        Screen('DrawTexture',window,...
                            obj.Init.txtrPointers{curEpoch}(curTxtrInd),...
                            [],obj.screenDim,[],0);
                        
                        % photodiode
                        Screen('FillRect', window, ...
                            pdContrast*obj.RGBscale, obj.pdDim); 
                        
                        % Flip the window at the framerate 
                        Screen('Flip', window, startTime + frameCount*ifi);   
                        pause(0) 
                        
                        % save current values of degShift and epoch
                        degShift(frameCount) = ...
                            obj.Init.epochDegShift{curEpoch}(curTxtrInd);
                        whichStim(frameCount) = 3; % moving is 3
                        whichEpoch(frameCount) = curEpoch;
                        
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
                obj.Out.whichStim = whichStim(1:(frameCount-1));
                obj.Out.whichEpoch = whichEpoch(1:(frameCount-1));
            else
                obj.Out.degShift = degShift;
                obj.Out.whichStim = whichStim;
                obj.Out.whichEpoch = whichEpoch;          
            end
        end 
        
        % Reconstructs the stimulus presented, for analysis
        function [rcStim, rcStimInd] = reconstructStim(obj, ...
                stimEpochStartTimes, lightStartTimes, darkStartTimes, ...
                order)
            display('Processing photodiode data for StaticAndMovingSquareGrating...');
            
            % lightStartTimes represent start of each new stimulus
            % length of flash (dark-light) defines which stimulus
            % call of pd transitions misses last light to dark transition
            if (length(darkStartTimes) < length(lightStartTimes))
                flashLengths = darkStartTimes - ...
                    lightStartTimes(1:length(darkStartTimes));
                endOnDark = 1;
            else
                endOnDark = 0;
            end
            % convert flashLengths to which stimulus presented
            stimPresented = round(flashLengths/obj.IFI);
            % 2 flashes corresponds to Out.whichStim code 1, 1 flash
            % corresponds to Out.whichStim code 2 - swap
            code1Ind = find(stimPresented==2);
            code2Ind = find(stimPresented==1);
            stimPresented(code1Ind) = 1;
            stimPresented(code2Ind) = 2;
            
            % if miss last stimulus transition, guess what stimPresented
            %  should be (+1 of previous), with looping
            if endOnDark
                stimPresented = [stimPresented; ...
                    mod(stimPresented(end),3) + 1];
            end
            
            % get what code thought it presented
            % to take care of 0 as last value in whichStim, whichEpoch
            if (obj.Out.whichStim(end) == 0)
                whichStimInd = [1; find(diff(obj.Out.whichStim(1:(end-1))))+1];
                whichStim = obj.Out.whichStim(whichStimInd);
            else
                whichStimInd = [1; find(diff(obj.Out.whichStim))+1];
                whichStim = obj.Out.whichStim(whichStimInd);   
            end
            if (obj.Out.whichEpoch(end) == 0)
                whichEpochInd = [1; find(diff(obj.Out.whichEpoch(1:(end-1))))+1];
                whichEpoch = obj.Out.whichEpoch(whichEpochInd);
            else
                whichEpochInd = [1; find(diff(obj.Out.whichEpoch))+1];
                whichEpoch = obj.Out.whichEpoch(whichEpochInd);
            end
            
            % reconstructed stimuli presented should match Out.whichStim
            % if whichStim does not equal stimPresented, we've missed a 
            %  stimulus transition or added a fake transition, don't 
            %  reconstruct in that case
            if (length(whichStim) ~= length(stimPresented))
                display('Cannot reconstruct stimulus.');
                rcStim = [];
                rcStimInd = [];
                return
            % two are same length (good!), check if values match
            else
                % if stimulus calls don't all match
                if ~(sum(whichStim==stimPresented) == length(whichStim))
                    display('PD signal and whichStim not the same. Using whichStim.');
                    % assume that single light frame was dropped, not that
                    %  entire stimulus was not presented; i.e. use
                    %  whichStim as correct
                    stimPresented = whichStim;                    
                end
            end
            
            % use stimPresented and whichEpoch to get transition times for
            % epochs and stimuli and which epoch was happening
            
            % stimulus transition times
            stimTransTimes = lightStartTimes;
            obj.Out.stimTransTimes = stimTransTimes;
            obj.Out.stimPresented = stimPresented;
            
            % epoch transition times
            epochTransTimes = lightStartTimes(stimPresented==1);
            obj.Out.epochTransTimes = epochTransTimes;
            obj.Out.epochsPresented = whichEpoch;
            
            % which epoch was being presented on each stimulus transition
            stimWhichEpoch = zeros(size(stimPresented));
            currEpochInd = 0;
            for i=1:length(stimPresented)
                % stimulus 1, gray, starts new epoch
                if (stimPresented(i) == 1)
                    currEpochInd = currEpochInd + 1;                    
                end
                    stimWhichEpoch(i) = whichEpoch(currEpochInd);
            end
            obj.Out.stimWhichEpoch = stimWhichEpoch;
            
            % transition times of individual stimulus types (not gray, as
            %  that's same as epoch transition)
            staticInd = find(stimPresented == 2);
            obj.Out.staticTransTimes = stimTransTimes(staticInd);
            obj.Out.staticWhichEpoch = stimWhichEpoch(staticInd);
            
            % moving
            movingInd = find(stimPresented == 3);
            obj.Out.staticTransTimes = stimTransTimes(movingInd);
            obj.Out.staticWhichEpoch = stimWhichEpoch(movingInd);
            
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
                rcStimInd(i) = whichEpoch(whichInd);
                % normalize to max of 1
                rcStim(i) = rcStimInd(i)/max(whichEpoch);
            end
             
        end
        
        % selects imaging frames to be used as F0 in computing dF/F during
        %  analysis
        % for this stimulus, gray period
        function [baselineFrameTimes, baselineSignals] = ...
                selectDFFbaseline(obj,imgFrameTimes, bksSignal, ...
                lightStartTimes, darkStartTimes)
            
            % initialize output variables
            baselineFrameTimes = [];
            baselineSignals = [];
            
            % get gray durations as array
            grayDur = cell2mat(obj.GrayDuration);
            
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
                endBaseTime = grayDur(whichEpoch); % end of gray
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