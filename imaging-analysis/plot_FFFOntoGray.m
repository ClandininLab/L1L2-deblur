% plot_FFFOntoGray.m
%
% Function to plot short flash off of gray responses. Plots each flash 
%  length's response separately. Plots all paired epochs onto 1 plot.
%  Expects 1D array (i.e. for 1 stimulus). Don't give this functioin
%  multiple columns for multiple stimuli
%
% INPUT:
%   rois - 1D array of structs, 1 for each ROI; must have field rats for
%       mean responses of that ROI, cell array of 1 for each epoch
%   inv - inverted? 1 for yes (negative goes up), 0 for no
%   yScale - y-axis limits
%   titleStr - genotype description to put in title
%   pairedEpochs - matrix indicating which epochs are paired (epochs in same
%       row are paired)
%
% OUTPUT:
%   creates plot as side effect, returns nothing
%
% Updated: 1/17/17 - address change in loadROIData for handling stimDat
% Updated: 5/18/17 HHY - change in compute_meanResp_moveAvg_FFFoG for
%  wrapping bins and adding t; change plotting to match
%
function plot_FFFOntoGray(rois,inv,yScale, titleStr, pairedEpochs)

    cm = colormap('lines');
    
%     flashDurs = cell2mat(rois(1).stimDat.obj.FlashDuration);
%     grayDur = rois(1).stimDat.obj.GrayDuration{1}; 
    
    % change in stimDat
    flashDurs = cell2mat(rois(1).stimDat.FlashDuration);
    grayDur = rois(1).stimDat.GrayDuration{1}; 
    
    for i = 1:size(pairedEpochs, 1) % number of flash lengths
        figure; 
        
        iDur = pairedEpochs(i);
        seqDur = flashDurs(iDur) + grayDur;

        for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
            rats = zeros(length(rois), length(rois(1).rats{pairedEpochs(i,j)}));
            flyID = zeros(1, length(rois));
            
            for n = 1:length(rois)
                rats(n, :) = rois(n).rats{pairedEpochs(i,j)};
                flyID(n) = rois(n).flyID;
            end    
            
            % time vector (use from 1st roi (all same))
            t = rois(1).t{pairedEpochs(i,j)};
                
        % get means and standard errors
            m = mean(rats,1);
            e = std(rats,[],1)./sqrt(size(rats,1));

%             plot_err_patch_v2(...
%                     (0:length(m))* rois(n).BIN_SHIFT,...
%                     [m m(1)],...
%                     [e e(1)], ...
%                     cm(j,:),(cm(j,:)+1)/2);
            plot_err_patch_v2(t,m,e,cm(j,:),(cm(j,:)+1)/2);
            hold on;
        end 
        % plotting parameters
        xScale = [0, seqDur];
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        xlim(xScale);
        line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line

        line([flashDurs(iDur) flashDurs(iDur)], yScale, 'color',[0 0 0],...
            'linestyle','--');
%         title([num2str(flashDurs(iDur)), ' sec flash']);
        legend('dark std error', 'dark', 'light stderror', 'light');
        if (inv)
            set(gca,'YDir','reverse');
        end
    
        titlename = [titleStr ', N cells = ' num2str(length(rois)), ...
            ', N flies = ', num2str(length(unique(flyID)))];
        title(titlename);
%         suplabel(titlename, 't');
        
        %         ', N flies = ' num2str(length(unique(flyID)))]);

    end
end  