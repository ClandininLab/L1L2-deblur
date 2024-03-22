% compute_SIMPLE_meanResp_err_moveAvg_NatStim_pool.m
%
% Function to compute simple mean response and error for naturalistic stimulus
%  for single ROI but multiple trials for that ROI. Works on 1 or more 
%  epochs (different stimulus sequences). All trials must have same number
%  of epochs
%  No resampling or moving average.
%
% INPUT
%   roi - struct(s) for a single ROI, which is one or more elements in the 
%       r x s ROI matrix (output of loadROIData)
%   seqDurs - length in seconds of each nat stim epoch
%   BIN_SHIFT - how much bins are separated in time - NOT USED CURRENTLY
%   BIN_WIDTH - how wide is a bin to average over - NOT USED CURRENTLY
%
% OUTPUT:
%   meanResp - mean response of this ROI to each epoch; cell array where
%       each element is different epoch
%   stdErr - std error associated with these responses; like meanResp, also
%       cell array
%
% Created: 2/16/22 - HHY
%
function [meanResp, stdErr, actualTimes, relativeTimes] = ...
    compute_simple_meanResp_err_moveAvg_NatStim_pool(roi, seqDurs, ...
    BIN_SHIFT, BIN_WIDTH)

    % number of epochs - assumes same number of epochs for trial
    nEpochs = length(roi(1).stimDat.WhichStim);
    
    % cells for saving 
    meanResp = {};
    stdErr = {};
    actualTimes = {}; %MMP
    relativeTimes = {}; %MMP

    % process each epoch type (diff. Nat Stim sequences)
    for i = 1:nEpochs  

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
            
            % save all frames after start time within this trial (MMP)
            for k = 1:length(thisEpochStartTimes)
                currStimStartTime = thisEpochStartTimes(k);
                currImFrameStartInd = find(roi(r).imFrameStartTimes-currStimStartTime > 0, 1);
                currImFrameInds = (currImFrameStartInd:currImFrameStartInd+framesPerEpoch-1);
                %this finds the first frame occuring after the start of the stim

                if length(roi(r).imFrameStartTimes)<currImFrameInds(end)
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
        end  
    end