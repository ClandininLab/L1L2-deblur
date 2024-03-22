% getRawSignalsLineBkgdSub.m
%
% Function that computes average fluorescence over time in all ROIs, with
%  and without background subtraction.
% Unlike getRawSignals.m, background subtraction here is line-by-line (Y
%  dimension). For each line, the average value in the background pixels
%  (as defined by the input mask) is subtracted from each ROI pixel in that 
%  line. This corrects for bleedthrough that isn't uniform across the whole
%  image
%
% INPUT:
%   in - struct of all of the imaging and stimulus values for a particular
%       time series
%   numLinesPool - number of lines to pool over for estimating background;
%       should divide evenly into number of lines (Y dim), otherwise, will 
%       use nearest integer
%   minPx - min background pixels per line(s) (specified by numLinesPool)
%       for estimating bkgd value 
%
% OUTPUT:
%   avSignal - average pixel intensity in each ROI, vector length k masks
%   bksSignal - average ROI intensity minus average background intensity
%
% CREATED: 8/3/22 - HHY
%
% UPDATED:
%   8/3/22 - HHY
%
function [avSignal, bksSignal] = getRawSignalsLineBkgdSub(in, ...
    numLinesPool, minPx)

    % number of masks
    nMasks = in.nMasks; 

    % get size of each frame
    numLines = size(in.alSeries,1);
    numPxPerLine = size(in.alSeries,2);

    % get number of lines to pool over, must divide evenly into numLines
    %  i.e. if it doesn't divide evenly, find nearest integer that does
    if (mod(numLines,numLinesPool))
        numLinesDiv = divisors(numLines); % divisors of numLines

        % get nearest divisor in value, set as numLinesPool
        [~, nearInd] = min(abs(numLinesDiv - numLinesPool));
        numLinesPool = numLinesDiv(nearInd);
    end

    % number of strips to analyze for each frame (total number of lines
    %  divided by numLinesPool)
    numStrips = numLines / numLinesPool;

    % check that minPx is 1 or more and less than or equal to numPxPerLine
    %  if not, set to 1
    if (minPx < 1) || (minPx > numPxPerLine)
        minPx = 1;
    end

    % preallocate arrays for ROI time series
    avSignal = zeros(nMasks, in.nFrames);
    bksSignal = zeros(nMasks, in.nFrames);
    sumSignal = zeros(nMasks, in.nFrames);
    bksSumSignal = zeros(nMasks, in.nFrames);
    sizeROIMasks = zeros(nMasks,1);
    masksInd = cell(nMasks,1);

    % get number of pixels in each ROI
    for k = 1:nMasks
        % number of pixels in each ROI
        sizeROIMasks(k) = sum(sum(in.roiMasks{k}));
        % indicies corresponding to pixels in ROI
        masksInd{k}= find(in.roiMasks{k}); 
    end

    % loop through all frames
    for i = 1:in.nFrames
        % this image frame (convert from uint16 to double)
        thisImg = double((in.alSeries(:,:,i)));

        % array for background value for each strip
        bkgdVals = zeros(numStrips,1);

        % loop through all strips
        for j = 1:numStrips
            % get indices into img rows of this strip
            thisStripStartInd = (j-1) * numLinesPool + 1;
            thisStripEndInd = thisStripStartInd + numLinesPool - 1;
            
            % logical mask for this strip
            thisStripMask = false(numLines, numPxPerLine);
            thisStripMask(thisStripStartInd:thisStripEndInd,:) = true;
            
            % this image frame, with all but this strip masked out
%             thisImgStrip = thisImg(thisStripMask);
            thisImgStrip = thisImg;
            thisImgStrip(~thisStripMask) = 0;

            % background
            % this image background region only
            thisImgStripBkgd = thisImgStrip(in.bkMask);

            % get number of background pixels in this strip
            numBkgdPx = sum(sum(thisStripMask & in.bkMask));

            % if number of background pixels in this strip is less than the
            %  minimum, use value in previous strip
            % if this is the first strip and the number of pixels is
            %  greater than 1, then use just these bkgd pixels
            % if there are no bkgd pixels for this first strip, don't
            %  background subtract for this strip
            if (numBkgdPx < minPx) % too few bkgd pixels
                if (j==1) % first strip
                    if (numBkgdPx < 1) % no background pixels
                        bkgdVals(j) = 0; % use 0
                    else % are background pixels, though less than min
                        % bkgd val is average across thse pixels
                        bkgdVals(j) = sum(thisImgStripBkgd) / numBkgdPx;
                    end
                else % not first strip
                    % bkgd val is val in previous strip
                    bkgdVals(j) = bkgdVals(j-1);
                end
            else % enough bkgd pixels
                bkgdVals(j) = sum(thisImgStripBkgd) / numBkgdPx;
            end

            % loop through all masks
            for k = 1:nMasks
                % ROI region within this strip
                maskImg = thisImgStrip(in.roiMasks{k}); 

                % number of pixels in this strip of ROI
                thisNumROIPx = sum(sum(thisStripMask & in.roiMasks{k}));

                % get sum of intensity in this part of ROI, cumulative over
                %  strips
                sumSignal(k,i) = sumSignal(k,i) + sum(maskImg);

                % sum of intensity in this part of ROI, with background
                %  subtraction subtract bkgd average * number of ROI pixels
                % cumulative over strips
                bksSumSignal(k,i) = bksSumSignal(k,i) + sum(maskImg) - ...
                    bkgdVals(j) * thisNumROIPx;
            end
        end
    end

    % get average intensity for each ROI
    % loop through all ROIs
    for k = 1:nMasks
        % average signal without background subtraction
        avSignal(k,:) = sumSignal(k,:) / sizeROIMasks(k);

        % average signal with background subtraction
        bksSignal(k,:) = bksSumSignal(k,:) / sizeROIMasks(k);
    end
end