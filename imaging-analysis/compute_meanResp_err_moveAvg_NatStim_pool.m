% compute_meanResp_err_moveAvg_NatStim_pool.m
%
% Function to compute mean response and error for naturalistic stimulus
%  for single ROI but multiple trials for that ROI. Works on 1 or more 
%  epochs (different stimlus sequences). All trials must have same number
%  of epochs
%  Computes mean response as moving average.
%
% INPUT
%   roi - struct(s) for a single ROI, which is one or more elements in the 
%       r x s ROI matrix (output of loadROIData)
%   seqDurs - length in seconds of each nat stim epoch
%   BIN_SHIFT - how much bins are separated in time
%   BIN_WIDTH - how wide is a bin to average over
%   FRAME_RATE - frame rate to interpolate to
%
% OUTPUT:
%   meanResp - mean response of this ROI to each epoch; cell array where
%       each element is different epoch
%   stdErr - std error associated with these responses; like meanResp, also
%       cell array
%
% Created: 2/16/22 - HHY
%
function [meanResp, stdErr] = ...
    compute_meanResp_err_moveAvg_NatStim_pool(roi, seqDurs, ...
    BIN_SHIFT, BIN_WIDTH, FRAME_RATE)

    % number of epochs - assumes same number of epochs for trial
    nEpochs = length(roi(1).stimDat.WhichStim);
    
    % cells for saving 
    meanResp = {};
    stdErr = {};

    % process each epoch type (diff. Nat Stim sequences)
    for i = 1:nEpochs  

        % pre-allocate for getting frame times relative to stimulus
        %  transition
        relTimes = [];
        dFF = [];

        % loop over all trials for this ROI (separate structs)
        for r = 1:length(roi)
        
            framesPerEpoch = ceil(seqDurs(i) * FRAME_RATE);
            
            % get info from stimulus reconstruction
            % epoch start times
            epochStartTimes = roi(r).pStimDat.rcStimInd;
            % contrast values for each epoch, as cell array
            epochIntVals = roi(r).pStimDat.rcStim;
            
            % split epoch start times by which epoch
            whichEpochs = roi(r).stimDat.epochsPresented;
            % clip epoch calls for those that were presented after aquisition
            %  stopped
            whichEpochs = whichEpochs(1:length(epochStartTimes));
            
            % start times just for this epoch type
            thisEpochStartTimes = epochStartTimes(whichEpochs==i);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save time when each imaging frame occured, relative to stimulus
            %  transition, and during which phase
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % calculate and save times relative to stimulus transition that
            %  imaging frame occured; save also dF/F
            for k = 1:length(roi(r).imFrameStartTimes)
    
                ft = roi(r).imFrameStartTimes(k); % current frame time
    
                % time difference between stimulus transitions and current
                %  frame time
                tDiff = ft - thisEpochStartTimes;
    
                % remove negative values - stim transition happened after frame
                tDiffValid = tDiff(tDiff>=0);
    
                % minimum value and index in time differences
                minTimeDiff = min(tDiffValid);
    
                % save relative time frame occured and dF/F value of that
                %  frame (for computing moving average)
                % ignore frames that happened in slightly long epochs
                if (minTimeDiff <= seqDurs(i)) 
                    relTimes = [relTimes minTimeDiff];
                    dFF = [dFF roi(r).dFF(k)];
                end               
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute moving average
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % pre-allocate
        slidingMean = zeros(1,framesPerEpoch);
        slidingStdErr = slidingMean;

        % number of bins that encroach around start and end 
        numEdgeBins = ceil((BIN_WIDTH/BIN_SHIFT - 1)/2);
        % number of bins in epoch
        numBins = framesPerEpoch;

        % compute sliding average
        for j = 1:numBins

            binCenter = (j-1)*BIN_SHIFT;
            binStart = binCenter - BIN_WIDTH/2;
            binEnd = binCenter + BIN_WIDTH/2;

            % at start, bin is smaller than usual
            if (j <= numEdgeBins) 
                binStart = 0;
            % at end, bin is smaller than usual    
            elseif (j > (numBins - numEdgeBins))
                binEnd = seqDurs(i);                
            end
            
%             valDFF = dFF((relTimes > binStart) & (relTimes <= binEnd));
            
            valInd = find((relTimes>binStart) .* (relTimes<=binEnd));
            valDFF = dFF(valInd);

            slidingMean(j) = mean(valDFF);
            slidingStdErr(j) = std(valDFF)/sqrt(length(valDFF));
        end
        
%         save output
        meanResp{i} = slidingMean;
        stdErr{i} = slidingStdErr; 
        
        % butterworth filter
%         d = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',24,'HalfPowerFrequency2',26, ...
%                'DesignMethod','butter','SampleRate',framesPerEpoch/seqDurs(i));
%         meanResp{i} = filtfilt(d,slidingMean);           
    end  
end