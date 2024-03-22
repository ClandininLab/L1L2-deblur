% filterROIs_SIMPLE_stimLockedPlots_NatStim_pool.m 
%
% Function to obtain stimulus-locked average responses and error per ROI
%  for presentations of Naturalistic Stimulus. Uses full 
%  contrast flash stimulus (search stimulus or full field) as way of 
%  identifying responding cells 
% Asks for user to identify responding/inverted/non-responding ROIs with
%  Y/I/N
% Works if there are multiple trials of the same stimulus (calls 
%  compute_meanResp_err_moveAvg_NatStim_pool.m)
%
% INPUT:
%   roiMat - 2D matrix of structs for each ROI, output of loadROIData
%   refColumn - which column of roiMat contains reference stimulus time
%       series (binary contrast or search stimulus); NOTE: this doesn't
%       really work. Assumes reference is in first column
%   FRAME_RATE - frame rate to bin average responses at, in Hz, or -1 to
%       skip resampling
%   binWidthMult - how many bins (at FRAME_RATE) to perform sliding average
%
% OUTPUT:
%   roiMatNew - roiMat input with additional fields rats and stdErr(s)
%       added, 1st column is ref stim, 2nd column is nat stim across all 
%       trials
%   iResp - indicies into roiMat of responding ROIs
%   iInv - indicies into roiMat of ROIs responding with opposite sign
%   framesPerLDCycle - number of frames in light flash + dark flash, for
%       reference stimulus
%
% Created: 2/16/22 - HHY
% 220607 MMP changed to skip resampling (using compute_simple_meanResp_...)

function [roiMatNew, iResp, iInv, framesPerLDCycle] = ...
    filterROIs_simple_stimLockedPlots_NatStim_pool(roiMat, refColumn,  ...
    FRAME_RATE, binWidthMult)
    
    % MMP get real frame rate
    if (FRAME_RATE == -1)
        FRAME_RATE = 1/median(diff(roiMat(1,1).imFrameStartTimes));
    end
    % initialization
    BIN_WIDTH = binWidthMult * (1/FRAME_RATE);
    BIN_SHIFT = 1/FRAME_RATE;
    P_VAL_THRESH = 0.01; % emperically, based on tests on fake data
    iResp = [];
    iInv = [];
    
    % timing for reference stimulus
%     refEpochDur = roiMat(1,1).stimDat.obj.Duration{1}*2;
    % change how stimDat is dealt with after loadROIData change
    refEpochDur = roiMat(1,refColumn).stimDat.Duration{1}*2;
    % timing for FFFoG stim is done in the for loop
    
    cm = colormap('lines');
    close all
    
    % number of epochs for naturalistic stimulus
    numEpochs = length(roiMat(1,2).stimDat.Duration);

    % only 2 columns, one for ref, one for nat stim
    roiMatNew = roiMat(:,1:2);
    
    % loop over all cells
    for r = 1:size(roiMat, 1)

        % get response to reference stimulus, in ref column

        % compute mean response, std err, indicies of
        %  significantly different regions of trace
        [refMeanResp, refStdErr, ~, ~, ~, ~, ~, framesPerLDCycle] = ...
        compute_meanResp_err_moveAvg_FFFBC(roiMat(r,refColumn), ...
        refEpochDur,BIN_SHIFT, BIN_WIDTH, FRAME_RATE, P_VAL_THRESH);
    
        % first col for ref
        roiMatNew(r,1).rats = refMeanResp;
        roiMatNew(r,1).stdErr = refStdErr;


        % get response for nat stim
        numCols = size(roiMat, 2);
        % get all columns that aren't refColumn
        natStimCols = 1:numCols;
        natStimCols(natStimCols == refColumn) = [];

        thisRoiNatStim = roiMat(r,natStimCols);

        % get length of epochs            
        epochDurations = cell2mat(thisRoiNatStim(1).stimDat.Duration);
        
        % compute mean and error - MMP changed to skip resampling
        [meanResp, stdErr] = ...
            compute_simple_meanResp_err_moveAvg_NatStim_pool(...
            thisRoiNatStim,epochDurations, BIN_SHIFT, BIN_WIDTH);
        
        % save info, 2nd col for nat stim
        roiMatNew(r,2).rats = meanResp;
        roiMatNew(r,2).stdErrs = stdErr;
        roiMatNew(r,2).BIN_SHIFT = BIN_SHIFT;


        % ask if reference stimulus is responding/not/inverted (but display
        % others)
        
        % display
        figure;

        % -----  plot reference stimulus ---- %
        subplot(size(roiMatNew,2), numEpochs, 1:numEpochs)
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
        for j = 2:size(roiMatNew,2)
            for k = 1:numEpochs
                subplot(size(roiMatNew,2), numEpochs,...
                    k+(j-1)*numEpochs)

                % plot response with error bar
                plot_err_patch_v2(...
                    (0:length(roiMatNew(r,j).rats{k}))...
                    * BIN_SHIFT,...
                    [roiMatNew(r,j).rats{k}...
                    roiMatNew(r,j).rats{k}(1)],...
                    [roiMatNew(r,j).stdErrs{k}...
                    roiMatNew(r,j).stdErrs{k}(1)], ...
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

%         close all
        
    end

end