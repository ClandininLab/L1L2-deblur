% plot_NatStimFFF.m
%
% Function to plot responses to naturalistic full field stimulus. One
%  figure per epoch (diff. seq). 2 subplots: neural response on top,
%  stimulus on bottom.
% Expects 1D array (i.e. for 1 stimulus). Don't give this function
%  multiple columns for multiple stimuli
%
% INPUT:
%   rois - 1D array of structs, 1 for each ROI; must have field rats for
%       mean responses of that ROI, cell array of 1 for each epoch
%   inv - inverted? 1 for yes (negative goes up), 0 for no
%   yScale - y-axis limits
%   titleStr - genotype description to put in title
%
% OUTPUT:
%   creates plot as side effect, returns nothing
%
% Created: 7/15/21 - HHY
%
function plot_NatStimFFF(rois,inv,yScale, titleStr)

    cm = colormap('lines');
    
    % number of epochs = number of separate figures
    numEpochs = length(rois(1).stimDat.WhichStim);
    
    % epoch durations
    epochDurations = cell2mat(rois(1).stimDat.Duration);
    
    % pull stimulus times and intensity values from reconstructed
    %  stimulus
    epochStimInt = rois(1).pStimDat.rcStim;
    whichEpochs = rois(1).stimDat.epochsPresented;
    
    for i = 1:numEpochs
        figure;
        
        
        % get means and standard errors
        % preallocate arrays
        rats = zeros(length(rois), length(rois(1).rats{i}));
        flyID = zeros(1, length(rois));
        
        % extract each cell's response and fly ID and put in array
        for n = 1:length(rois)
            rats(n,:) = rois(n).rats{i};
            flyID(n) = rois(n).flyID;
        end
            
        % get means and standard errors
        m = mean(rats,1);
        e = std(rats,[],1)./sqrt(size(rats,1));
        
        % plot neuronal response
        subplot(2,1,1); % 2 horizontal plots, top plot
        plot_err_patch_v2(...
                (0:length(m))* rois(n).BIN_SHIFT,...
                [m m(1)],...
                [e e(1)], ...
                cm(i,:),(cm(i,:)+1)/2);
            hold on;
        % plotting parameters
        xScale = [0, epochDurations(i)];
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        xlim(xScale);
        ylim(yScale);
        line([0 epochDurations(i)],[0 0],'color',[0 0 0]); % x-axis line
        if (inv)
            set(gca,'YDir','reverse');
        end
    
        titlename = [titleStr ', N cells = ' num2str(length(rois)), ...
            ', N flies = ', num2str(length(unique(flyID)))];
        title(titlename);
        
        % plot stimulus
        % find intensity sequence corresponding to this epoch
        ind = find(whichEpochs(2:end) == i) + 1;
        stimSeq = epochStimInt{ind};
        stimT = linspace(0, epochDurations(i), length(stimSeq));
        
        subplot(2,1,2); % 2 horizontal plots, bottom plot
        plot(stimT, stimSeq);
        
        % plotting parameters
        xScale = [0, epochDurations(i)];
        xlabel('time (sec)');
        ylabel('Intensity');
        xlim(xScale);
        ylim([0 1]);
        

    end
end  