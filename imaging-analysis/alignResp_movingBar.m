% alignResp_movingBar.m
%
% Function to align responses to moving bar, to allow averaging across
%  cells.
% User selects peak in x and y. Shift to center point.
%

function roiMat = alignResp_movingBar(roiMat,slctPtFrac)

    % loop through all ROIs
    for r = 1:length(roiMat)
        % loop through all directions
        for i = 1:length(roiMat(r).rats)
            figure;
            plot(roiMat(r).rats{i});
            title('Select peak');
            [peakInd,~] = ginput(1); % get user selection 
            
            centerInd = slctPtFrac*length(roiMat(r).rats{i});
            posShift = round(centerInd - peakInd);
            
            roiMat(r).shiftRats{i} = circshift(roiMat(r).rats{i},...
                [0,posShift]);
            
            close;
        end
        
    end

end