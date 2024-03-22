% getRawSignals.m
% 
% Function that computes average fluorescence over time in all the ROIs.
%  Outputs both average fluorescence with and without background
%  subtraction.
%  The functions drawROIs and selectBackground must have been run first to
%  define the ROIs and background region
%
% INPUT:
%   in - struct of all of the imaging and stimulus values for a particular
%       time series
%
% OUTPUT:
%   avSignal - average pixel intensity in each ROI, vector length k masks
%   bksSignal - average ROI intensity minus average background intensity
%

function [avSignal, bksSignal] = getRawSignals(in)

    nMasks = in.nMasks; 

    % preallocate arrays for ROI time series
    avSignal = zeros(nMasks, in.nFrames);
    bksSignal = zeros(nMasks, in.nFrames);
    sizeROIMasks = zeros(nMasks,1);
    
    for k = 1:nMasks
        % number of pixels in each ROI
        sizeROIMasks(k) = sum(sum(in.roiMasks{k}));
        % indicies corresponding to pixels in ROI
        masksInd{k}= find(in.roiMasks{k}); 
    end

    % as above, but for background, not ROIs
    sizeBkMask = sum(sum(in.bkMask));
    % bkMaskInd = find(in.bkMask);

    % compute signal in each ROI
    for i = 1:in.nFrames
        A = double((in.alSeries(:,:,i))); % individual frame image
        bkMaskImg = A(in.bkMask); % background region

        for k = 1:nMasks
            maskImg = A(masksInd{k}); % ROI region
            % average intensity in ROI (sum over space, divide by ROI size)
            avSignal(k,i) = sum(maskImg)./sizeROIMasks(k);
            % background subtraction (by average background intensity)
            bksSignal(k,i) = avSignal(k,i) - sum((bkMaskImg))./sizeBkMask;
        end
    end

end 