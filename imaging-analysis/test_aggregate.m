% test_aggregate.m
%
% Function to obtain stimulus-locked average response and error per ROI 
%  and which ROIs are responding, for full field flash stimulus
% Average response and error per ROI computed with moving average, with
%  averaging window of size BIN_WIDTH and shift between averaging windows
%  of BIN_SHIFT
% To select responding ROIs, displays trace and computes p-values comparing
%  binned time points in dark flash vs. light flash;in green, mean by 
%  binning; in blue, mean by sliding windown average
%
% INPUT:
%    in - 1D array of structs (indexed by ROI), output of createROIMatrix
%    cycleDur - length (duration in seconds) of light + dark flash
%    FRAME_RATE - in Hz, for bin shift
%    binWidthMult - width of bin as this multiplier of 1/FRAME_RATE. If you
%       did not upsample or interpolate you imaging frames, binWidthMult=1
%
% OUTPUT:
%   out - struct that outputs stimulus-locked average response and error,
%     arrayed by ROI, plus which ROIs are responding or inverted, plus 
%     some metadata 
%   respInd - indicies of ROIs that are responding
%   invInd - indicies of ROIs that are responding with inverted sign
%   framesPerCycle - number of frames in light + dark flash
%   BIN_WIDTH - width of bin, in seconds
%   P_VAL_THRESH - p-value threshold for determining whether bins in
%       corresponding times during light and dark flash are significantly
%       different (currently, not being used)
%
% 160629 update: the code that uses the results of the ttest are commented
% out because it encounters an error during the ttest if we interpolate the
% responses 
%
function [out, respInd, invInd, framesPerCycle, BIN_WIDTH, P_VAL_THRESH] = ...
    test_aggregate(in, cycleDur, FRAME_RATE, binWidthMult)
    
    BIN_WIDTH = binWidthMult * (1/FRAME_RATE);
    BIN_SHIFT = 1/FRAME_RATE;
    
    P_VAL_THRESH = 0.01; % emperically, based on tests on fake data
    
    % number of frames per stimulus L+D epoch
%     framesPerEpoch = ceil(epochDur * FRAME_RATE); % valid?

     % pre-allocate output arrays
    out = struct([]);
    respInd = [];
    invInd = [];
    
    for r = 1:length(in) % loops through all ROIs
        % pass in this ROI's data to compute mean response and std error
        % across light/dark trials
        [meanResp, stdErr, ~, meanDark, meanLight, stdErrDark, stdErrLight, ...
            framesPerCycle]...
            = compute_meanResp_err_moveAvg_FFFBC(in(r), cycleDur,...
            BIN_SHIFT, BIN_WIDTH, FRAME_RATE, P_VAL_THRESH);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display moving average mean trace and significant pVals
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xScale = [0, cycleDur * 2];
        figure; hold on;
    
% %         debug
%            length(meanResp)
%            length(meanDark)
%            length((0:length(meanResp)*2)*BIN_SHIFT)
%             length( [meanDark meanLight meanDark meanLight meanDark(1)])

        % plot in green binned average
        plot_err_patch_v2((0:length(meanResp)*2)*BIN_SHIFT,...
            [meanDark meanLight meanDark meanLight meanDark(1)],...
            [stdErrDark stdErrLight stdErrDark stdErrLight stdErrDark(1)],...
            [0 1 0],[0.5 1 0.5]);
        
        % debug
        %  size(meanResp)
        %  size(stdErr)
        
        % plot in blue moving average
        plot_err_patch_v2((0:length(meanResp)*2)*BIN_SHIFT,...
            [meanResp; meanResp; meanResp(1)], ...
            [stdErr; stdErr; stdErr(1)], [0 0 1],...
            [0.5 0.5 1]);
        xlabel('time (sec)');ylabel('response (dF/F)');
        yScale = ylim;
        xlim(xScale);
        line([0 cycleDur*2],[0 0],'color',[0 0 0]);
        for k = 1:4
            line([(k-1)*cycleDur/2 (k-1)*cycleDur/2],yScale,'color',[0 0 0],'linestyle','--');
        end
        patch([0 0 cycleDur/2 cycleDur/2],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)],...
            [0 0 0]);
        patch([cycleDur/2 cycleDur/2 cycleDur cycleDur],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)]...
            ,[1 1 1]);
        patch([cycleDur cycleDur cycleDur/2*3 cycleDur/2*3],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)],...
            [0 0 0]);
        patch([cycleDur/2*3 cycleDur/2*3 cycleDur*2 cycleDur*2],...
            [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)]...
            ,[1 1 1]);
        set(gca,'xTick',0:(cycleDur/2):(cycleDur*2));
%         set(gca,'YDir','reverse')
        title(sprintf('Cell %d, from fly ID %d of %s', r, ...
            in(r).flyID, in(r).seriesID));

        % place a red * for all bins where p-value < P_VAL_THRESH
        
%         % convert significant bins to times
%         sigPIndDup = [sigPInd, sigPInd + length(darkVals), ...
%             sigPInd + 2 * length(darkVals), sigPInd + 3 * length(darkVals)];
%         sigPIndTimes = (sigPIndDup-1) / FRAME_RATE + (1/FRAME_RATE/2); 
%         
%         plot(sigPIndTimes,ones(1,length(sigPIndTimes))*yScale(2) * 0.75,...
%             'r*');
        
        
        % get user input on whether cell is responding or not
        userInput = input('Responding? [Y/I/N]: ','s');
        
        if (strcmpi(userInput, 'y')) % yes, responding
            respInd = [respInd r];
        elseif (strcmpi(userInput, 'i')) % inverted
            invInd = [invInd r];
        end % anything else, treat as no
        
        close all % close current figure after input
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save in output arrays
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        out(r).rats = meanResp;
        out(r).stdErr = stdErr;
        out(r).seriesID = in(r).seriesID;
        out(r).cellNum = in(r).roiMask;
        out(r).flyID = in(r).flyID;
    end
    
end