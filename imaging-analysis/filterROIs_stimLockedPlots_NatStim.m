% filterROIs_stimLockedPlots_NatStim.m 
%
% Function to obtain stimulus-locked average responses and error per ROI
%  for presentations of Naturalistic Stimulus. Uses full 
%  contrast flash stimulus (search stimulus or full field) as way of 
%  identifying responding cells 
% Asks for user to identify responding/inverted/non-responding ROIs with
%  Y/I/N
%
% INPUT:
%   roiMat - 2D matrix of structs for each ROI, output of loadROIData
%   refColumn - which column of roiMat contains reference stimulus time
%       series (binary contrast or search stimulus)
%   FRAME_RATE - frame rate to bin average responses at, in Hz
%   binWidthMult - how many bins (at FRAME_RATE) to perform sliding average
%
% OUTPUT:
%   roiMat - roiMat input with additional fields rats and stdErr(s) added
%   iResp - indicies into roiMat of responding ROIs
%   iInv - indicies into roiMat of ROIs responding with opposite sign
%   framesPerLDCycle - number of frames in light flash + dark flash, for
%       reference stimulus
%
% Created: 7/14/21 - HHY
%

function [roiMat, iResp, iInv, framesPerLDCycle] = ...
    filterROIs_stimLockedPlots_NatStim(roiMat, refColumn,  ...
    FRAME_RATE, binWidthMult)
    
    % initialization
    BIN_WIDTH = binWidthMult * (1/FRAME_RATE);
    BIN_SHIFT = 1/FRAME_RATE;
    P_VAL_THRESH = 0.01; % emperically, based on tests on fake data
    iResp = [];
    iInv = [];
    
    % timing for reference stimulus
%     refEpochDur = roiMat(1,1).stimDat.obj.Duration{1}*2;
    % change how stimDat is dealt with after loadROIData change
    refEpochDur = roiMat(1,1).stimDat.Duration{1}*2;
    % timing for FFFoG stim is done in the for loop
    
    cm = colormap('lines');
    close all
    
    % number of epochs for naturalistic stimulus
    numEpochs = length(roiMat(1,2).stimDat.Duration);
    
    % loop over all cells
    for r = 1:size(roiMat, 1)
        for s = 1:size(roiMat,2)
            % if ROI is from ref stim, just plot response using values
            % saved in the data matrix 
            if (s == refColumn)
                % compute mean response, std err, indicies of
                %  significantly different regions of trace
                [refMeanResp, refStdErr, ~, ~, ~, ~, ~, framesPerLDCycle] = ...
                compute_meanResp_err_moveAvg_FFFBC(roiMat(r,s), refEpochDur,...
                BIN_SHIFT, BIN_WIDTH, FRAME_RATE, P_VAL_THRESH);
            
                roiMat(r,s).rats = refMeanResp;
                roiMat(r,s).stdErr = refStdErr;

                % if not reference stimulus
            else
                % get length of epochs            
                epochDurations = cell2mat(roiMat(r,s).stimDat.Duration);
                
                % compute mean and error
                [meanResp, stdErr] = ...
                    compute_meanResp_err_moveAvg_NatStim(roiMat(r,s),...
                    epochDurations, BIN_SHIFT, BIN_WIDTH, FRAME_RATE);
                
                % save info
                roiMat(r,s).rats = meanResp;
                roiMat(r,s).stdErrs = stdErr;
                roiMat(r,s).BIN_SHIFT = BIN_SHIFT;
%                     in(r,s).sigPInd = sigPInd;
            end 
        end 
        % ask if reference stimulus is responding/not/inverted (but display
        % others)
        
        % display
        figure;

        % -----  plot reference stimulus ---- %
        subplot(size(roiMat,2), numEpochs, 1:numEpochs)
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
        
        % ----- plot naturalistic stimuli ---- % 
        for j = 2:size(roiMat,2)
            for k = 1:numEpochs
                subplot(size(roiMat,2), numEpochs,...
                    k+(j-1)*numEpochs)

                % plot response with error bar
                plot_err_patch_v2(...
                    (0:length(roiMat(r,j).rats{k}))...
                    * BIN_SHIFT,...
                    [roiMat(r,j).rats{k}...
                    roiMat(r,j).rats{k}(1)],...
                    [roiMat(r,j).stdErrs{k}...
                    roiMat(r,j).stdErrs{k}(1)], ...
                    cm(k,:),(cm(k,:)+1)/2);

                seqDur = epochDurations(k);

                xScale = [0, seqDur];
                xlabel('time (sec)');
                ylabel('response (dF/F)');
                yScale = ylim;
                xlim(xScale);
                line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line

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