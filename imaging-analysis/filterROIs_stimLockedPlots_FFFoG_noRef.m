% filterROIs_stimLockedPlots_FFFoG_noRef.m
%
% Function to obtain stimulus-locked average responses and error per ROI
%  for presentations of different length flashes contrast levels off of
%  gray. Does not expect search stimulus or binary contrast reference
%  stimulus. Instead, works on any number of matched flash off gray
%  stimuli, or single one.
%
% INPUT:
%   roiMat - ROI data structs, output of loadROIData or createROIMatrix
%   pairedEpochs - 2D matrix indicating which epochs of flash off gray
%       stimulus are paired with each other (e.g. light and dark of same 
%       flash duration)
%   FRAME_RATE - frame rate to bin average responses at, in Hz
%   binWidthMult - how many bins (at FRAME_RATE) to perform sliding average
%   
% OUTPUT:
%   roiMat - roiMat input with additional fields rats and stdErr(s) added
%   iResp - indicies into roiMat of responding ROIs
%   iInv - indicies into roiMat of ROIs responding with opposite sign
%
% Updated: 4/21/17 - MX - to include code for inverting the responses in
%   the plots if inv = 1 (for ASAP imaging)

function [roiMat, iResp, iInv] = ...
    filterROIs_stimLockedPlots_FFFoG_noRef(roiMat, pairedEpochs, ...
    FRAME_RATE, inv, yScale, binWidthMult)
    
    % initialization
    BIN_WIDTH = binWidthMult * (1/FRAME_RATE);
    BIN_SHIFT = 1/FRAME_RATE;
    P_VAL_THRESH = 0.01; % emperically, based on tests on fake data
    iResp = [];
    iInv = [];

    % timing for FFFoG stim is done in the for loop
    
    cm = colormap('lines');
    close all
    
    % for placing stimulus patch appropriately with inverted or not
    if (inv)
    	patchY = yScale(1);
    else
        patchY = yScale(2);
    end
    
    % loop over all cells
    for r = 1:size(roiMat, 1)
        % reject right away if moving 
        for s = 1:size(roiMat,2)
            % if ROI is from ref stim, just plot response using values
            % saved in the data matrix 
%             if (s == refColumn)
%                 % compute mean response, std err, indicies of
%                 %  significantly different regions of trace
%                 [refMeanResp, refStdErr] = ...
%                 compute_meanResp_err_moveAvg_FFFBC(roiMat(r,s), refEpochDur,...
%                 BIN_SHIFT, BIN_WIDTH, FRAME_RATE, P_VAL_THRESH);
%             
%                 roiMat(r,s).rats = refMeanResp;
%                 roiMat(r,s).stdErr = refStdErr;
% 
%                 % if not reference stimulus
%             elseif (s == fffgColumn)
                % get length of flash epochs and gray epochs 
                flashDurations = cell2mat(roiMat(r,s).stimDat.FlashDuration);
                grayDuration = roiMat(r,s).stimDat.GrayDuration{1}; 
                
                % compute mean and error
                [meanResp, stdErr,t] = ...
                    compute_meanResp_err_moveAvg_FFFoG(roiMat(r,s),...
                    flashDurations + grayDuration, ...
                    BIN_SHIFT, BIN_WIDTH, FRAME_RATE);
                
                % save info
                roiMat(r,s).rats = meanResp;
                roiMat(r,s).stdErrs = stdErr;
                roiMat(r,s).t = t;
                roiMat(r,s).BIN_SHIFT = BIN_SHIFT;
%                     in(r,s).sigPInd = sigPInd;
%             end 
        end 
        % ask if reference stimulus is responding/not/inverted (but display
        % others)
        
        % display
        figure;

        % ----- plot full field flash onto gray stimuli ---- % 
        % for each light/dark pairing 
        for j = 1:size(roiMat,2)
            for k = 1:size(pairedEpochs,1)
                subplot(size(roiMat,2),size(pairedEpochs,1),...
                    k+(j-1)*size(pairedEpochs,1))

                for l = 1:size(pairedEpochs,2)
                    plot_err_patch_v2(...
                        roiMat(r,j).t{pairedEpochs(k,l)},...
                        roiMat(r,j).rats{pairedEpochs(k,l)},...
                        roiMat(r,j).stdErrs{pairedEpochs(k,l)},...
                        cm(l,:),(cm(l,:)+1)/2);
                    hold on;
                end
                iDur = pairedEpochs(k);
                seqDur = flashDurations(iDur) + grayDuration;

                xScale = [0, seqDur];
                xlabel('time (sec)');
                ylabel('response (dF/F)');
                yScale = ylim;
                xlim(xScale);
                line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line

                line([flashDurations(iDur) flashDurations(iDur)], yScale, 'color',[0 0 0],...
                    'linestyle','--');
                
                if (inv) % reverse axes on inverted
                     set(gca,'YDir','reverse');
                end
            end
        end
        
        % ask if responding
        userInput = input('Responding? [Y/I/N]: ','s');
        
        if (strcmpi(userInput, 'y')) % yes, responding
            iResp = [iResp r];
        elseif (strcmpi(userInput, 'i')) % inverted
            iInv = [iInv r];
        end % anything else, treat as no

        close all
        
    end

end