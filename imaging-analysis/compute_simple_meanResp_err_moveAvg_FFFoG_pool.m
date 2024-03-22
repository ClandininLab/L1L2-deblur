% compute_SIMPLE_meanResp_err_moveAvg_FFFoG_pool.m
%
% Function to compute simple mean response and error for different length 
%   flashes off of gray stimulus for single ROI. No resampling or moving
%   average.
%
% INPUT
%   roi - struct(s) for a single ROI, which is one or more elements in the 
%       r x s ROI matrix (output of loadROIData)
%   seqDurs - length in seconds of each flash+gray epoch (vector with 1
%       value for each epoch)
%   BIN_SHIFT - how much bins are separated in time - NOT USED CURRENTLY
%   BIN_WIDTH - how wide is a bin to average over - NOT USED CURRENTLY
%   FRAME_RATE - frame rate to interpolate to - NOT USED CURRENTLY
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
    compute_simple_meanResp_err_moveAvg_FFFoG_pool(roi, seqDurs, ...
    BIN_SHIFT, BIN_WIDTH, FRAME_RATE)    

    % number of epochs - assumes same number of epochs for trial
    nFlashEpochs = length(roi(1).stimDat.FlashContrast);
    
    % cells for saving 
    meanResp = {};
    stdErr = {};
    t = {}; % time stamps on meanResp and stdErr
    actualTimes = {}; %MMP
    relativeTimes = {}; %MMP

    
    % process each type of flash-to-gray sequence
    for i = 1:nFlashEpochs  
        % pre-allocate for getting frame times relative to stimulus
        %  transition
        relTimes = [];
        actTimes = [];
        dFF = [];

        % loop over all trials for this ROI (separate structs)
        for r = 1:length(roi)
        
            % set the epoch duration, rounding down (MMP)
            real_frame_rate = 1/median(diff(roi(r).imFrameStartTimes));
            framesPerEpoch = floor(seqDurs(i) * real_frame_rate);
    
            % gray epoch to flash epoch times 
            rcStimInd = roi(r).pStimDat.rcStimInd;
            g2fTimes = roi(r).pStimDat.stimEpochStartTimes(rcStimInd == i);
            % (start times just for this epoch type)

            % save all frames after start time within this trial (MMP)
            for k = 1:length(g2fTimes)
                currStimStartTime = g2fTimes(k);
                currImFrameStartInd = find(roi(r).imFrameStartTimes-currStimStartTime > 0, 1);
                currImFrameInds = (currImFrameStartInd:currImFrameStartInd+framesPerEpoch-1);
                %this finds the first frame occurring after the start of the stim
                
                if isempty(currImFrameInds)
                    break
                elseif length(roi(r).imFrameStartTimes)<currImFrameInds(end)
                    break
                end                
                actTimes = [actTimes roi(r).imFrameStartTimes(currImFrameInds)];
                relTimes = [relTimes roi(r).imFrameStartTimes(currImFrameInds)-currStimStartTime];
                dFF = [dFF roi(r).dFF(currImFrameInds)];
            end
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute straight average WITHOUT MOVING AVERAGE (MMP)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            meanResp{i} = (mean(dFF,2))';
            stdErr{i} = (std(dFF,0,2)/sqrt(length(dFF)))';
            actualTimes{i} = actTimes;
            relativeTimes{i} = relTimes;
            t{i} = mean(relativeTimes{i},2);
        end
    end