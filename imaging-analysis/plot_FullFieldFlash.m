% plot_FullFieldFlash.m
%
% Function to plot full field flash responses. Plots light+dark flash
%  responses 2X, so all transitions are visible.
% Works on any binary contrast flashing stimulus (e.g. binary contrast full
%  field flash and search stimulus)
%
% INPUT:
%   in - array of ROI structs, has field rats for mean response of each
%       cell
%   inv - 1 for inverted (e.g. V-Ind), 0 for not inverted axes (i.e.
%       inverted is negative going up)
%   yScale - y-axis limits
%   ifi - interframe interval (in sec)
%   stimTime - length of single flash (in sec)
%   titleStr - genotype description to put in title
%   indivCells - 1 for plotting individual cells behind mean, 0 for no
%
% OUTPUT:
%   none but generates plot as side effect
%
function plot_FullFieldFlash(in, inv, yScale, ifi, ...
    stimTime, titleStr, indivCells)
    rats = [in.rats]; % rows are values, columns are ROIs
    nROIs = size(rats,2); % number of ROIs
    framesPerCycle = size(rats,1);
    
    xScale = [0, stimTime*4];
%     framesPerEpoch = ceil(stimTime*2*(1/ifi));

    % key plotting values
    m = mean(rats,2);
    e = std(rats, 0, 2)./sqrt(nROIs);
    t = [0:(framesPerCycle*2)]'*ifi;

    % debug
%     length(t)
%     size(rats)
%     size([rats(:,1); rats(:,1); rats(1,1)])
    
    % for placing stimulus patch appropriately with inverted or not
    if (inv)
    	patchY = yScale(1);
    else
        patchY = yScale(2);
    end
    
    figure; hold on;
    
    % Plot traces of individual cells
    if (indivCells)
        for i=1:length(in)
%                 length([rats(:,i); rats(:,i); rats(1,i)]) % debug
        plot(t,[rats(:,i); rats(:,i); rats(1,i)], 'color', ...
            [0.7 0.7 0.7]);
        end
    end
    
    % Plot the mean response and the std error
    h1 = plot_err_patch_v2(t,[m; m; m(1)],...
        [e; e; e(1)],[0 0 1],[0.5 0.5 1]);
    
    xlabel('time (sec)');ylabel('response (dF/F)');
    ylim(yScale);
    xlim(xScale);
    line([0 stimTime*4],[0 0],'color',[0 0 0]);
    for k = 1:4
        line([(k-1)*stimTime (k-1)*stimTime],yScale,'color',[0 0 0],'linestyle','--');
    end
    patch([0 0 stimTime stimTime],...
        [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
        [0 0 0]);
    patch([stimTime stimTime stimTime*2 stimTime*2],...
        [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
        ,[1 1 1]);
    patch([stimTime*2 stimTime*2 stimTime*3 stimTime*3],...
        [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
        [0 0 0]);
    patch([stimTime*3 stimTime*3 stimTime*4 stimTime*4],...
        [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
        ,[1 1 1]);
    set(gca,'xTick',0:stimTime:stimTime*4);
    if (inv) % reverse axes on inverted
        set(gca,'YDir','reverse');
    end
    title([titleStr ', N cells = ' num2str(nROIs) ...
        ', flies = ' num2str(length(unique([in.flyID]')))]);
    niceaxes;
end