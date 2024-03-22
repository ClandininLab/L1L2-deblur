% aggregate_fff_matchedROIs.m 
%
% Function to obtain stimulus-locked average responses and error per ROI
%  for different presentations of full field flash. 
% Resamples to specified frame rate. Allows boxcar averaging.
% Uses full contrast flash stimulus as way of identifying responding
%  cells (either binary contrast or search stimulus).
% Asks users whether each ROI is responding, inverted, or not-responding.
%  Displays plot, asks users for Y/N/I.
% Adaptation of filterROIs_stimLockedPlots_FFFoG and test_aggregate
%
% INPUT:
%   roiMat - 2D array of structs, rows are ROIs, columns are different time
%       series; output of loadROIData and then sortROIDataMat function
%   refColumn - numerical index indicating which column is reference column 
%   FRAME_RATE - frame rate to bin average responses at
%   binWidthMult - how many bins (at FRAME_RATE) to perform sliding average
%
% OUTPUT:
%   roiMat - input roiMat, with fields rats and stdErr (for ROI's average
%       response and the standard error of that response) added
%   iResp - indices into rows of roiMat, which ROIs are responding (user
%       flagged)
%   iInv - indices into rows of roiMat, which ROIs are responding with
%       inverted sign from expected (user flagged) 
%   BIN_WIDTH - width of each bin of ROI's average response, determined
%       from FRAME_RATE and binWidthMult
%
% Updated: 3/13/17 - HHY - created to analyze matched fff data
% Updated: 3/17/17 - HHY - added comments
%

function [roiMat, iResp, iInv, BIN_WIDTH] = ...
    aggregate_fff_matchedROIs(roiMat, refColumn, FRAME_RATE, binWidthMult)
    
    % initialization
    BIN_WIDTH = binWidthMult * (1/FRAME_RATE);
    BIN_SHIFT = 1/FRAME_RATE;
    P_VAL_THRESH = 0.01; % emperically, based on tests on fake data
    iResp = [];
    iInv = [];
    
    refEpochDur = roiMat(1,1).stimDat.Duration{1}*2;
    
    % loop over all cells
    for r = 1:size(roiMat, 1)
        % reject right away if moving 
        for s = 1:size(roiMat,2)
            if (s == refColumn)
                % compute mean response, std err, indicies of
                %  significantly different regions of trace
                [refMeanResp, refStdErr, ~, ~, ~, ~, ~, ~] = ...
                compute_meanResp_err_moveAvg_FFFBC(roiMat(r,s), ...
                refEpochDur,BIN_SHIFT, BIN_WIDTH, FRAME_RATE, P_VAL_THRESH);

                roiMat(r,s).rats = refMeanResp;
                roiMat(r,s).stdErr = refStdErr;
            else
                % compute mean response, std err, indicies of
                %  significantly different regions of trace
                [meanResp, stdErr, ~, ~, ~, ~, ~, ~] = ...
                compute_meanResp_err_moveAvg_FFFBC(roiMat(r,s), ...
                roiMat(r,s).stimDat.Duration{1}*2,...
                BIN_SHIFT, BIN_WIDTH, FRAME_RATE, P_VAL_THRESH);

                roiMat(r,s).rats = meanResp;
                roiMat(r,s).stdErr = stdErr;        
            end
        end 
        % ask if reference stimulus is responding/not/inverted (but display
        % others)
        
        % display
        figure;

        % -----  plot reference stimulus ---- %
        subplot(2,size(roiMat,2)-1,1:(size(roiMat,2)-1))
        title(sprintf('Cell %d',r));
        xScale = [0, refEpochDur * 2];
        plot_err_patch_v2((0:length(refMeanResp)*2)*BIN_SHIFT,...
            [refMeanResp; refMeanResp; refMeanResp(1)], ...
            [refStdErr; refStdErr; refStdErr(1)], [0 0 1],[0.5 0.5 1]);
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        yScale = ylim;
        xlim(xScale);
        line([0 refEpochDur*2],[0 0],'color',[0 0 0]);
        for k = 1:4
            line([(k-1)*refEpochDur/2 (k-1)*refEpochDur/2],yScale,...
                'color',[0 0 0],'linestyle','--');
        end
        patch([0 0 refEpochDur/2 refEpochDur/2],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)],...
            [0 0 0]);
        patch([refEpochDur/2 refEpochDur/2 refEpochDur refEpochDur],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)]...
            ,[1 1 1]);
        patch([refEpochDur refEpochDur refEpochDur/2*3 refEpochDur/2*3],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)],...
            [0 0 0]);
        patch([refEpochDur/2*3 refEpochDur/2*3 refEpochDur*2 refEpochDur*2],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)]...
            ,[1 1 1]);
        set(gca,'xTick',0:(refEpochDur/2):(refEpochDur*2));
        
        % plot red * for all bins where p-value < P_VAL_THRESH
%         sigPIndTimes = (sigPInd-1) / FRAME_RATE + (1/FRAME_RATE/2);
%         plot(sigPIndTimes,ones(1,length(sigPIndTimes))*yScale(1) * 0.75,...
%             'r*');
        
        % ----- plot non-reference ---- % 
        % for each column
        for k = 2:size(roiMat,2)
            subplot(2,size(roiMat,2)-1,k+size(roiMat,2)-2)

            plot_err_patch_v2((0:length(roiMat(r,k).rats)*2)*BIN_SHIFT,...
                [roiMat(r,k).rats; roiMat(r,k).rats; roiMat(r,k).rats(1)], ...
                [roiMat(r,k).stdErr; roiMat(r,k).stdErr; ...
                roiMat(r,k).stdErr(1)], [0 0 1],[0.5 0.5 1]);
            epochDur = roiMat(r,k).stimDat.Duration{1}*2;
            xlabel('time (sec)');
            ylabel('response (dF/F)');
            yScale = ylim;
            xlim(xScale);
            line([0 epochDur*2],[0 0],'color',[0 0 0]);
            for j = 1:4
                line([(j-1)*epochDur/2 (j-1)*epochDur/2],yScale,...
                    'color',[0 0 0],'linestyle','--');
            end
            patch([0 0 epochDur/2 epochDur/2],...
                [(yScale(2)*0.85) (yScale(2)*0.9) ...
                (yScale(2)*0.9) (yScale(2)*0.85)],...
                [0 0 0]);
            patch([epochDur/2 epochDur/2 epochDur epochDur],...
                [(yScale(2)*0.85) (yScale(2)*0.9) ...
                (yScale(2)*0.9) (yScale(2)*0.85)]...
                ,[1 1 1]);
            patch([epochDur epochDur epochDur/2*3 epochDur/2*3],...
                [(yScale(2)*0.85) (yScale(2)*0.9) ...
                (yScale(2)*0.9) (yScale(2)*0.85)],...
                [0 0 0]);
            patch([epochDur/2*3 epochDur/2*3 epochDur*2 epochDur*2],...
                [(yScale(2)*0.85) (yScale(2)*0.9) ...
                (yScale(2)*0.9) (yScale(2)*0.85)]...
                ,[1 1 1]);
            set(gca,'xTick',0:(epochDur/2):(epochDur*2));

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