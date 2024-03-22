% computePolarInd_staticMovingGrating.m
%
% Function that finds peak response during static and moving periods of
% static and moving gratings stimulus for each epoch. For making polar
% plots. Adds to ROI data struct.
%
% 

function roiMat = computePolarInd_staticMovingGrating(roiMat, frameRate)

    % loop through all cells
    for r = 1:length(roiMat)
        % get static duration, gray duration
        staticDuration = cell2mat(roiMat(r).stimDat.StaticDuration);
        grayDuration = cell2mat(roiMat(r).stimDat.GrayDuration);
        
        % preallocate
        peakStatic = zeros(size(staticDuration));
        peakMoving = zeros(size(staticDuration));
        
        % loop through all epochs
        for i = 1:length(staticDuration)
            % find static peak 
            staticStartInd = ceil(grayDuration(i)*frameRate);
            staticEndInd = floor((grayDuration(i)+staticDuration(i))*...
                frameRate);
            peakStatic(i) = max(...
                roiMat(r).rats{i}(staticStartInd:staticEndInd));
            
            % find moving peak
            movingStartInd = ceil((grayDuration(i)+staticDuration(i))*...
                frameRate);
            peakMoving(i) = max(roiMat(r).rats{i}(movingStartInd:end));
            
        end
        
        % save
        roiMat(r).peakStatic = peakStatic;
        roiMat(r).peakMoving = peakMoving;
    end

end