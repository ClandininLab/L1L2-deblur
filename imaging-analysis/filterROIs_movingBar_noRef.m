% filterROIs_movingBar_noRef.m
%
% Function to obtain stimulus-locked average responses and error per ROI
%  for presentations of moving bar stimulus. No reference stimulus
%
% INPUT:
%   roiMat - ROI data structs, output of loadROIData or createROIMatrix
%   FRAME_RATE - frame rate to bin average responses at, in Hz
%   binWidthMult - how many bins (at FRAME_RATE) to perform sliding average
%   
% OUTPUT:
%   roiMat - roiMat input with additional fields rats and stdErr(s) added
%   iResp - indicies into roiMat of responding ROIs
%   iInv - indicies into roiMat of ROIs responding with opposite sign
%   layer - layer ROI belongs to
%  
% Updated: 4/25/17 - still WIP, HHY

function [roiMat, iResp, iInv, layer] = ...
    filterROIs_movingBar_noRef(roiMat, FRAME_RATE, binWidthMult)
    
    % initialization
    BIN_WIDTH = binWidthMult * (1/FRAME_RATE);
    BIN_SHIFT = 1/FRAME_RATE;
    P_VAL_THRESH = 0.01; % emperically, based on tests on fake data
    iResp = [];
    iInv = [];
    layer = [];

    % timing for FFFoG stim is done in the for loop
    
    cm = colormap('lines');
    close all
    
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
                % get length of epochs
                epochDuration = cell2mat(roiMat(r,s).stimDat.EpochDuration);
                
                % compute mean and error
                [meanResp, stdErr] = ...
                    compute_meanResp_err_moveAvg_movingBar(roiMat(r,s),...
                    epochDuration, BIN_SHIFT, BIN_WIDTH, FRAME_RATE);
                
                % save info
                roiMat(r,s).rats = meanResp;
                roiMat(r,s).stdErrs = stdErr;
                roiMat(r,s).BIN_SHIFT = BIN_SHIFT;
%                     in(r,s).sigPInd = sigPInd;
%             end 
        end 
        % ask if reference stimulus is responding/not/inverted (but display
        % others)
        
        % display
        figure;

        % -----  plot reference stimulus ---- %
%         subplot(2,size(pairedEpochs,1),1:size(pairedEpochs,1))
%         title(sprintf('Cell %d',r));
%         xScale = [0, refEpochDur * 2];
%         plot_err_patch_v2((0:length(refMeanResp)*2)*BIN_SHIFT,...
%             [refMeanResp; refMeanResp; refMeanResp(1)], ...
%             [refStdErr; refStdErr; refStdErr(1)], [0 0 1],[0.5 0.5 1]);
%         xlabel('time (sec)');
%         ylabel('response (dF/F)');
%         yScale = ylim;
%         xlim(xScale);
%         line([0 refEpochDur*2],[0 0],'color',[0 0 0]);
%         for k = 1:4
%             line([(k-1)*refEpochDur/2 (k-1)*refEpochDur/2],yScale,...
%                 'color',[0 0 0],'linestyle','--');
%         end
%         patch([0 0 refEpochDur/2 refEpochDur/2],...
%             [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)],...
%             [0 0 0]);
%         patch([refEpochDur/2 refEpochDur/2 refEpochDur refEpochDur],...
%             [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)]...
%             ,[1 1 1]);
%         patch([refEpochDur refEpochDur refEpochDur/2*3 refEpochDur/2*3],...
%             [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)],...
%             [0 0 0]);
%         patch([refEpochDur/2*3 refEpochDur/2*3 refEpochDur*2 refEpochDur*2],...
%             [(yScale(2)*0.85) (yScale(2)*0.9) (yScale(2)*0.9) (yScale(2)*0.85)]...
%             ,[1 1 1]);
%         set(gca,'xTick',0:(refEpochDur/2):(refEpochDur*2));
        
        % plot red * for all bins where p-value < P_VAL_THRESH
%         sigPIndTimes = (sigPInd-1) / FRAME_RATE + (1/FRAME_RATE/2);
%         plot(sigPIndTimes,ones(1,length(sigPIndTimes))*yScale(1) * 0.75,...
%             'r*');
        
        % ----- plot moving bar responses ---- % 
        % for each epoch
        numEpochs = length(roiMat(r).rats);
        for j = 1:numEpochs
            subplot(numEpochs,1,j);
            plot_err_patch_v2(...
                (0:length(roiMat(r).rats{j})) * BIN_SHIFT,...
                [roiMat(r).rats{j} roiMat(r).rats{j}(1)],...
                [roiMat(r).stdErrs{j} roiMat(r).stdErrs{j}(1)],...
                cm(1,:),(cm(1,:)+1)/2);
            hold on;
            
            seqDur = epochDuration(j);
            
            xScale = [0, seqDur];
            xlabel('time (sec)');
            ylabel('response (dF/F)');
            yScale = ylim;
            xlim(xScale);
            if (j == 1)
                title(['Cell ' num2str(r)]);
            end
            
            line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line
        end
        
        set(gcf, 'units','normalized','outerposition',[0.2 0 0.4 1]);

        
        % ask if responding
        userInput = input('Responding? [Y/I/N]: ','s');
        
        if (strcmpi(userInput, 'y')) % yes, responding
            iResp = [iResp r];
            
            % if responding, ask user for which layer cell is in
            load([char(roiMat(r).seriesID) '_pData.mat']);
            
            cm = colormap('lines');
            close

            figure;
            
            curColor = cm(1,:);
            curMask = cat(3,curColor(1).*pDat.roiMasks{roiMat(r).roiMask},...
                curColor(2).*pDat.roiMasks{roiMat(r).roiMask},...
                curColor(3).*pDat.roiMasks{roiMat(r).roiMask});
           
            imshow(pDat.avgIm,[]);
            hold on;
            h = imshow(curMask);
            set(h,'AlphaData',0.5); % make color mask slightly transparent
            
            userInput2 = input('Which layer? [1/2/3/4]: ');
            layer = [layer userInput2];
            
        elseif (strcmpi(userInput, 'i')) % inverted
            iInv = [iInv r];
        end % anything else, treat as no

        close all
        
    end

end