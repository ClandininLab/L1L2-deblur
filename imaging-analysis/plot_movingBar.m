% plot_movingBar.m
%
% To plot mean response to moving bar, 1 plot for each layer, 1 subplot for
%  each direction.

function plot_movingBar(respDataMat, iLayer, yScale, plotTitle, ...
    epochOrder, frameRate, epochDuration)
    % loop over each layer
    cm = colormap('lines');
    close all
    
    for i=1:4
        figure;
        % cells in this layer
        thisRespDataMat = respDataMat(iLayer == i);
        
        if ~isempty(thisRespDataMat)
            % loop over directions
            for k = 1:length(thisRespDataMat(1).shiftRats)
                shiftRats = zeros(length(thisRespDataMat),...
                    length(thisRespDataMat(1).shiftRats{k}));
                flyID = zeros(length(thisRespDataMat),1);
                % loop over all cells, save response, flyID
                for j=1:length(thisRespDataMat)
                    shiftRats(j,:) = thisRespDataMat(j).shiftRats{k};
                    flyID(j) = thisRespDataMat(j).flyID;
                end

                % compute mean, std err
                m = mean(shiftRats,1);
                e = std(shiftRats,[],1)/sqrt(length(thisRespDataMat));
                
                % plot
                subplot(length(thisRespDataMat(1).shiftRats),1,epochOrder(k));
                plot_err_patch_v2((0:length(m)) * (1/frameRate),[m m(1)],...
                    [e e(1)], cm(i,:),(cm(i,:)+1)/2);
                hold on;
                
                seqDur = epochDuration(k);
            
                xScale = [0, seqDur];
                xlabel('time (sec)');
                ylabel('response (dF/F)');
                ylim(yScale);
                xlim(xScale);
                if (k == 1)
                    title([plotTitle ' layer ' num2str(i), ...
                        ' n = ' num2str(length(thisRespDataMat)) '(', ...
                        num2str(length(unique(flyID))) ')']);
                end

                line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line
            end
            set(gcf, 'units','normalized','outerposition',[0.2 0 0.4 1]);
        end
        
    end
end