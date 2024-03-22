% plot_gratings_polar.m
%
% function to make polar plots from static and moving gratings stimulus, of
% OS and DS

function plot_gratings_polar(respDataMat, iLayer, yScale, plotTitle)

    cm = colormap('lines');
    close all
    
    % loop over each layer
    for i=1:4
        % cells in this layer
        thisRespDataMat = respDataMat(iLayer == i);
        
        if ~isempty(thisRespDataMat)
            peakStatic = zeros(length(thisRespDataMat),...
                length(thisRespDataMat(1).peakStatic));
            peakMoving = zeros(length(thisRespDataMat),...
                length(thisRespDataMat(1).peakMoving));
            flyID = zeros(length(thisRespDataMat),1);
            % loop over all cells, save response, flyID
            for j=1:length(thisRespDataMat)
                peakStatic(j,:) = thisRespDataMat(j).peakStatic;
                peakMoving(j,:) = thisRespDataMat(j).peakMoving;
                flyID(j) = thisRespDataMat(j).flyID;
            end
            
            % compute mean
            mStatic = mean(peakStatic,1);
            eStatic = std(peakStatic,[],1)/sqrt(length(thisRespDataMat));
            mMoving = mean(peakMoving,1);
            eMoving = std(peakMoving,[],1)/sqrt(length(thisRespDataMat));
            
            % plot
            thetaShift = 2*pi/length(thisRespDataMat(1).peakStatic);
            theta = 0:thetaShift:(2*pi);
            
            figure;
            polarwitherrorbar(theta,[mStatic mStatic(1)],...
                [eStatic eStatic(1)]);
            title([plotTitle ' static layer ' num2str(i) ' n = ',...
                num2str(length(thisRespDataMat)) '(' ,...
                num2str(length(unique(flyID))) ')']);
            
            figure;
            polarwitherrorbar(theta,[mMoving mMoving(1)],...
                [eMoving eMoving(1)]);
            title([plotTitle ' moving layer ' num2str(i) ' n = ',...
                num2str(length(thisRespDataMat)) '(' ,...
                num2str(length(unique(flyID))) ')']);
            
        end
    end
end