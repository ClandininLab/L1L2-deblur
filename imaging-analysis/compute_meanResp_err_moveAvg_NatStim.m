% compute_meanResp_err_moveAvg_NatStim.m
%
% Function to compute mean response and error for naturalistic stimulus
%  for single ROI. Works on 1 or more epochs (different stimlus sequences)
%  Computes mean response as moving average.
%
% INPUT
%   roi - struct for a single ROI, which is one element in the 
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
% Created: 7/14/21 - HHY
%
function [meanResp, stdErr] = ...
    compute_meanResp_err_moveAvg_NatStim(roi, seqDurs, ...
    BIN_SHIFT, BIN_WIDTH, FRAME_RATE)    
    
    % number of epochs
    nEpochs = length(roi.stimDat.WhichStim);
    
    % cells for saving 
    meanResp = {};
    stdErr = {};
    
    % process each epoch type (diff. Nat Stim sequences)
    for i = 1:nEpochs  
        
        framesPerEpoch = ceil(seqDurs(i) * FRAME_RATE);
        
        % get info from stimulus reconstruction
        % epoch start times
        epochStartTimes = roi.pStimDat.rcStimInd;
        % contrast values for each epoch, as cell array
        epochIntVals = roi.pStimDat.rcStim;
        
        % split epoch start times by which epoch
        whichEpochs = roi.stimDat.epochsPresented;
        % clip epoch calls for those that were presented after aquisition
        %  stopped
        whichEpochs = whichEpochs(1:length(epochStartTimes));
        
        % start times just for this epoch type
        thisEpochStartTimes = epochStartTimes(whichEpochs==i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Bin data for computing p-values
        % Save time when each imaging frame occured, relative to stimulus
        %  transition, and during which phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % pre-allocate for getting frame times relative to stimulus
        %  transition
        relTimes = [];
        dFF = [];
        
        % pre-allocate cell array for binning imaging frames
        binnedVals = cell(framesPerEpoch,1);
        
        % place each imaging frame in appropriate bin
        % calculate and save times relative to stimulus transition that
        %  imaging frame occured; save also dF/F
        for k = 1:length(roi.imFrameStartTimes)

            ft = roi.imFrameStartTimes(k); % current frame time

            % time difference between stimulus transitions and current
            %  frame time
            tDiff = ft - thisEpochStartTimes;

            % remove negative values - stim transition happened after frame
            tDiffValid = tDiff(tDiff>=0);

            % minimum value and index in time differences
            minTimeDiff = min(tDiffValid);

            % put imaging frame into appropriate bin
            whichBin = ceil(minTimeDiff * FRAME_RATE);
            if (whichBin == 0) % corrects for rare case where index=0
                whichBin = 1;
            end

            % only take imaging frames within duration of this epoch (and
            %  not belonging to another epoch)
            if (whichBin <= framesPerEpoch)
                binnedVals{whichBin} = [binnedVals{whichBin} roi.dFF(k)];
            end
            % save relative time frame occured and dF/F value of that
            %  frame (for computing moving average)
            % ignore frames that happened in slightly long epochs
            if (minTimeDiff <= seqDurs(i)) 
                relTimes = [relTimes minTimeDiff];
                dFF = [dFF roi.dFF(k)];
            end               
        end
        
        % compute mean and standard error by binning (not sliding)
        % currently not used, but could swap out in future
        
%         % pre-allocate arrays
%         meanBinning = zeros(1,framesPerEpoch);
%         stdErrBinning = meanBinning;
% 
%         for l = 1:(framesPerEpoch)
%             meanBinning(l) = mean(binnedVals{l});
%             stdErrBinning(l) = std(binnedVals{l})/sqrt(length(binnedVals{l}));
%         end
%         
%         % deal with NaNs (if any bins have no values)
%         [meanBinning, stdErrBinning] = correctShortFlashBinNaN(meanBinning,...
%             stdErrBinning);
        
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

        % output without sliding average
%         meanResp{i} = meanBinning;
%         stdErr{i} = stdErrBinning;         
        
%         binnedVals{i} = vals;
    end
    
%     % compute p-values
%     for i=1:size(pairedEpochs,1)
%         vals1 = binnedVals{pairedEpochs(i,1)};
%         vals2 = binnedVals{pairedEpochs(i,2)};
%         
%         pVals = ones(1,length(vals1));
%         for j=1:length(vals1) % num of bins
%             [h,p] = ttest2(vals1{j},vals2{j});
%             if(~isempty(p))
%                 pVals(j) = p;
%             end
%         end
%         
%         % save p-values
%         sigPInd = find(pVals < P_VAL_THRESH);
%         sigPIndDup{pairedEpochs(i,1)} = sigPInd;
%         sigPIndDup{pairedEpochs(i,2)} = sigPInd;  
%     end     
end