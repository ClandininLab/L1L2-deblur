% compute_meanResp_err_moveAvg_FFFBC.m
%
% Computes mean response and std error across trials for a single ROI
%  for FullFieldFlash_BinaryContrast or SearchStimulusFlash stimuli
% Resamples to specified frame rate by binning individual samples based on
%  when the frame was captured relative to a stimulus transition
% 
% INPUT:
%   roi - struct of info for single ROI, of the format created by
%       loadROIData and createROIMatrix
%   cycleDur - duration of light flash + dark flash
%   BIN_SHIFT - time between 2 adjacent bins in average
%   BIN_WIDTH - width of each bin (i.e. based on BIN_SHIFT and BIN_WIDTH,
%       they can be overlapping - boxcar averaging)
%   FRAME_RATE - frame rate average response is sampled to = 1/BIN_SHIFT
%   P_VAL_THRESH - p-value threshold for determining whether bins in
%       corresponding times during light and dark flash are significantly
%       different (currently, not being used)
%
% OUTPUT:
%   slidingMean - mean response of ROI, moving average
%   slidingStdErr- standard error on response of ROI, for computation of
%       mean by moving average
%   sigPInd - which indicies are where responses during light and dark are
%       significantly different
%   meanDark - mean response of ROI during dark period, just binned 
%       response, not moving average
%   meanLight - mean response of ROI during light period, just binned
%       response, not moving average
%   stdErrDark - standard error on meanDark
%   stdErrLight - standard error on meanLight
%   framesPerLDCycle - number of frames in light flash + dark flash; if
%       cycleDur isn't evenly divided by FRAME_RATE, this is rounded so it
%       is
%
function [slidingMean, slidingStdErr, sigPInd, meanDark, meanLight, ...
    stdErrDark, stdErrLight, framesPerLDCycle] = ...
    compute_meanResp_err_moveAvg_FFFBC(roi, ...
    cycleDur, BIN_SHIFT, BIN_WIDTH, FRAME_RATE, P_VAL_THRESH)

    framesPerLDCycle = ceil(cycleDur * FRAME_RATE);
    if mod(framesPerLDCycle, 2) > 0
        framesPerLDCycle = framesPerLDCycle - 1;
    end 
    framesPerEpoch = framesPerLDCycle/2;

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bin data for computing p-values
    % Save time when each imaging frame occured, relative to stimulus
    %  transition, and during which phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % find stimulus transition times
    % dark to light
    dtlTimes = roi.pStimDat.lightStartTimes;
    % light to dark
    ltdTimes = roi.pStimDat.darkStartTimes;

    % do not count imaging frames before first dark to light transition
    startImgFrame = find((roi.imFrameStartTimes - dtlTimes(1))>0,...
        1,'first');

    % pre-allocate for getting frame times relative to stimulus
    %  transition
    lightRelTimes = [];
    darkRelTimes = [];
    lightDFF = [];
    darkDFF = [];

    % pre-allocate cell arrays for binning imaging frames
    darkVals = cell(framesPerEpoch, 1);
    lightVals = darkVals;

    % place each imaging frame in appropriate bin
    % calculate and save times relative to stimulus transition that
    %  imaging frame occured; save also dF/F
    for k = startImgFrame:length(roi.imFrameStartTimes)
        % current frame time
        ft = roi.imFrameStartTimes(k); 
        % time difference between stimulus transitions and current
        %  frame time
        tDiff = ft - roi.pStimDat.stimEpochStartTimes;
        % remove negative values - stim transition happened after frame
        tDiffValid = tDiff(tDiff >= 0);
        % minimum value and index in time differences
        [minTimeDiff, ind]= min(tDiffValid);

        % put the imaging frame in the correct bin in the stimulus
        % epoch based on it's time relative to the stimulus transition
        % each bin is about the width of an imaging frame
        whichBin = ceil(minTimeDiff / (1/FRAME_RATE));
        if (whichBin == 0) % corrects for rare case where index=0
            whichBin = 1;
        end

        % stimulus contrast at that timestamp was light, put imaging 
        % frame in the light bin
        if roi.pStimDat.rcStim(ind) == 1 
            % when stimulus presented slightly longer than expected,
            %  ignore that img frame
            if (whichBin <= framesPerEpoch)
                lightVals{whichBin} = [lightVals{whichBin}; roi.dFF(k)];
            end
            % if time difference is less than or equal to length of
            % light epoch, save the dFF value
            if (minTimeDiff <= cycleDur/2) 
                lightRelTimes = [lightRelTimes; minTimeDiff];
                lightDFF = [lightDFF; roi.dFF(k)];
            end 
        else 
            % when stimulus presented slightly longer than expected,
            %  ignore that img frame
            if (whichBin <= framesPerEpoch)
                darkVals{whichBin} = [darkVals{whichBin}; roi.dFF(k)];
            end

            % save relative time frame occured and dF/F value of that
            %  frame (for computing moving average)
            % ignore frames that happened in slightly long flashes
            if (minTimeDiff <= cycleDur/2)
                darkRelTimes = [darkRelTimes; minTimeDiff];
                darkDFF = [darkDFF; roi.dFF(k)];
            end                
        end        
    end


    % compute mean and standard error by binning (not sliding)
    % pre-allocate arrays 
    meanDark = zeros(1, framesPerEpoch);
    meanLight = meanDark;
    stdErrDark = meanDark;
    stdErrLight = meanDark;

    for l = 1:(framesPerEpoch)
        meanDark(l) = mean(darkVals{l});
        stdErrDark(l) = std(darkVals{l})/sqrt(length(darkVals{l}));
        meanLight(l) = mean(lightVals{l});
        stdErrLight(l) = std(lightVals{l})/sqrt(length(lightVals{l}));
    end

    % deal with NaNs (if any bins have no values)
    [meanDark, stdErrDark] = correctFFFBinNaN(meanDark, meanLight,...
        stdErrDark);
    [meanLight, stdErrLight] = correctFFFBinNaN(meanLight, meanDark,...
        stdErrLight);

    sigPInd = [];
%         % compute p-values
%         pVals = zeros(1,length(darkVals));
%         for j=1:length(darkVals) % num of bins
%             %Debug
%             darkVals = darkVals{j} 
%             lightVals = lightVals{j}
%             
%             [h,p] = ttest2(darkVals{j},lightVals{j});
%             pVals(j) = p;
%         end
%         sigPInd = find(pVals < P_VAL_THRESH);
%         

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute moving average
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % pre-allocate
    totalBins = framesPerLDCycle;
    slidingMean = zeros(totalBins, 1);
    slidingStdErr = slidingMean;

    % number of bins that encroach around start and end 
    numEdgeBins = ceil((BIN_WIDTH/BIN_SHIFT - 1)/2);
    % number of bins in each flash
    numBinsFlash = totalBins/2;

    % compute sliding average, both dark and light; start with dark
    for j = 1:totalBins

        binCenter = (j-1)*BIN_SHIFT;
        binStart = binCenter - BIN_WIDTH/2;
        binEnd = binCenter + BIN_WIDTH/2;

        % wrap at start
        if (j <= numEdgeBins) 
            binStart = binStart + cycleDur/2;

            % obtain appropriate values from both light and dark data
            darkValDFF = darkDFF(darkRelTimes <= binEnd);
            lightValDFF = lightDFF(lightRelTimes > binStart);
            valDFF = [darkValDFF; lightValDFF];

        % at dark->light edge
        elseif ((j <= numBinsFlash) && ... 
                (j > (numBinsFlash - numEdgeBins)))
            binEnd = binEnd - cycleDur/2;

            % obtain appropriate values from both light and dark data
            darkValDFF = darkDFF(darkRelTimes > binStart);
            lightValDFF = lightDFF(lightRelTimes <= binEnd);
            valDFF = [darkValDFF; lightValDFF];         

        % at dark->light edge, corner case of right edge (includes 0
        %  in light times)
        elseif (j == (numBinsFlash - numEdgeBins))
            % obtain appropriate values from both light and dark data
            valInd = find((darkRelTimes > binStart) .* ...
                (darkRelTimes <= binEnd));
            darkValDFF = darkDFF(valInd);
            lightValDFF = lightDFF(lightRelTimes == 0);
            valDFF = [darkValDFF; lightValDFF];

        % wrap at start of light epoch
        elseif ((j > numBinsFlash) && (j <= numBinsFlash + numEdgeBins))
            binEnd = binEnd - (cycleDur/2);
            % binStart - no modification b/c - and + epochDur/2

            darkValDFF = darkDFF(darkRelTimes > binStart);
            lightValDFF = lightDFF(lightRelTimes <= binEnd);
            valDFF = [darkValDFF; lightValDFF];

        % at light->dark edge
        elseif (j > (numBinsFlash * 2 - numEdgeBins))
            binStart = binStart - (cycleDur/2);
            binEnd = binEnd - cycleDur;

            darkValDFF = darkDFF(darkRelTimes <= binEnd);
            lightValDFF = lightDFF(lightRelTimes > binStart);
            valDFF = [darkValDFF; lightValDFF];

        % at light->dark edge, corner case of right edge (includes 0 in
        %  dark times)
        elseif (j == (numBinsFlash * 2 - numEdgeBins))
            binStart = binStart - (cycleDur/2);
            binEnd = binEnd - (cycleDur/2);

            darkValDFF = darkDFF(darkRelTimes == 0);
            valInd = find((lightRelTimes > binStart) .* ... 
                (lightRelTimes <= binEnd));
            lightValDFF = lightDFF(valInd);
            valDFF = [darkValDFF; lightValDFF];

        % regular, during dark flash
        elseif (j <= numBinsFlash)
            valInd = find((darkRelTimes>binStart) .*...
                (darkRelTimes<=binEnd));
            valDFF = darkDFF(valInd);
        % all else, during light flash
        else
            binStart = binStart - (cycleDur/2);
            binEnd = binEnd - (cycleDur/2);
            valInd = find((lightRelTimes>binStart) .*...
                (lightRelTimes<=binEnd));
            valDFF = lightDFF(valInd);            
        end

        slidingMean(j) = mean(valDFF);
        slidingStdErr(j) = std(valDFF)/sqrt(length(valDFF));
        
    end       
       
end 