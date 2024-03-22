% compute_meanResp_err_moveAvg_FFFoG_pool.m
%
% Function to compute mean response and error for different length flashes
%   off of gray stimulus for single ROI. Computes mean response as moving
%   average.
%
% INPUT
%   roi - struct(s) for a single ROI, which is one or more elements in the 
%       r x s ROI matrix (output of loadROIData)
%   seqDurs - length in seconds of each flash+gray epoch (vector with 1
%       value for each epoch)
%   BIN_SHIFT - how much bins are separated in time
%   BIN_WIDTH - how wide is a bin to average over
%   FRAME_RATE - frame rate to interpolate to
%
% OUTPUT:
%   meanResp - mean response of this ROI to each epoch; cell array where
%       each element is different epoch
%   stdErr - std error associated with these responses; like meanResp, also
%       cell array
%   t - time stamps for bin center for meanResp and stdErr
%
% CREATED: 3/29/22 - HHY
%
function [meanResp, stdErr, t] = ...
    compute_meanResp_err_moveAvg_FFFoG_pool(roi, seqDurs, ...
    BIN_SHIFT, BIN_WIDTH, FRAME_RATE)    

    % number of epochs - assumes same number of epochs for trial
    nFlashEpochs = length(roi(1).stimDat.FlashContrast);
    
    % cells for saving 
    meanResp = {};
    stdErr = {};
    t = {}; % time stamps on meanResp and stdErr


    
    % process each type of flash-to-gray sequence
    for i = 1:nFlashEpochs  
        % pre-allocate for getting frame times relative to stimulus
        %  transition
        relTimes = [];
        dFF = [];

        % loop over all trials for this ROI (separate structs)
        for r = 1:length(roi)
        
            % set the epoch duration depending on whether its a flash or gray
    %         if i > nFlashEpochs % if gray epoch
    %             epochDur = roi.stimDat.obj.GrayDuration{1};
    %         else 
    %             epochDur = roi.stimDat.obj.FlashDuration{i};
    %         end 
            
    %         framesPerEpoch = ceil(seqDurs(i) * FRAME_RATE);
            % round down
            framesPerEpoch = floor(seqDurs(i) * FRAME_RATE);
    
            % gray epoch to flash epoch times 
            rcStimInd = roi(r).pStimDat.rcStimInd;
            g2fTimes = roi(r).pStimDat.stimEpochStartTimes(rcStimInd == i);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Bin data for computing p-values
            % Save time when each imaging frame occured, relative to stimulus
            %  transition, and during which phase
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    
            
            % pre-allocate cell array for binning imaging frames
            binnedVals = cell(framesPerEpoch,1);
            
            % start allocating frames into bins after first gray to flash
            %  transition
    %         startImgFrame = find((roi.imFrameStartTimes - g2fTimes(1))>=0,1,'first');
            
            % place each imaging frame in appropriate bin
            % calculate and save times relative to stimulus transition that
            %  imaging frame occured; save also dF/F
            for k = 1:length(roi(r).imFrameStartTimes)
    
                ft = roi(r).imFrameStartTimes(k); % current frame time
    
                % time difference between stimulus transitions and current
                %  frame time
                tDiff = ft - g2fTimes;
    
    %             % remove negative values - stim transition happened after frame
    %             tDiffValid = tDiff(tDiff>=0);
    % 
    %             % minimum value and index in time differences
    %             minTimeDiff = min(tDiffValid);
    
                % remove values where stim transition happened too far after
                %  frame (-1*BIN_WIDTH later)
                tDiffValid = tDiff(tDiff>=(-1*BIN_WIDTH));
                
                % minimum time difference relative to stimulus transition
                [~, indMin] = min(abs(tDiffValid));
                minTimeDiff = tDiffValid(indMin);
    
                % put imaging frame into appropriate bin
                whichBin = ceil(minTimeDiff * FRAME_RATE);
    %             if (whichBin == 0) % corrects for rare case where index=0
    %                 whichBin = 1;
    %             end
    
                % when stimulus presented slightly longer than expected,
                %  ignore that img frame; bin only values that happen after
                %  stimulus transition (whichBin > 0)
                if ((whichBin <= framesPerEpoch) & (whichBin > 0))
                    binnedVals{whichBin} = [binnedVals{whichBin} roi(r).dFF(k)];
                end
                % save relative time frame occured and dF/F value of that
                %  frame (for computing moving average)
                % ignore frames that happened in slightly long flashes
                if (minTimeDiff <= seqDurs(i))
                    relTimes = [relTimes minTimeDiff];
                    dFF = [dFF roi(r).dFF(k)];
                end               
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
        % length of array is framesPerEpoch + 1, to account for wrapping
        slidingMean = zeros(1,framesPerEpoch+1);
        slidingStdErr = slidingMean;
        movAvgTimes = slidingMean;

        % number of bins in epoch
        numBins = framesPerEpoch+1;

        % compute sliding average
        for j = 1:numBins
            binCenter = (j-1)*BIN_SHIFT - BIN_WIDTH/2;
            binStart = binCenter - BIN_WIDTH/2;
            binEnd = binCenter + BIN_WIDTH/2;

%             % at start, bin is smaller than usual
%             if (j <= numEdgeBins) 
%                 binStart = 0;
%             % at end, bin is smaller than usual    
%             elseif (j > (numBins - numEdgeBins))
%                 binEnd = seqDurs(i);                
%             end
            
%             valDFF = dFF((relTimes > binStart) & (relTimes <= binEnd));
            
            valInd = find((relTimes>binStart) .* (relTimes<=binEnd));
            valDFF = dFF(valInd);

            slidingMean(j) = mean(valDFF);
            slidingStdErr(j) = std(valDFF)/sqrt(length(valDFF));
            movAvgTimes(j) = binEnd;
        end
        
%         save output
        meanResp{i} = slidingMean;
        stdErr{i} = slidingStdErr; 
        t{i} = movAvgTimes;
        
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