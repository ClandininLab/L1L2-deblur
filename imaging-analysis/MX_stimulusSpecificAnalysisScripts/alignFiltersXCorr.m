% alignFiltersXCorr.m
% 
% Aligns filters in the time axis when filters are erroneously time-shifted 
% due to frame-dropping/adding in the LightCrafter. 
% This operation requires the past-future filters and won't work on the
% regular (past-present) filters.
%
% INPUTS: 
%   pfFilters - past/future linear filters for all cells 
%   iRef - index of the reference filte in pfFilters 
%   tau - half-length of the pfFilters in seconds
%   stimIFI - stimulus interframe interval after upsampling 
%   startCut - seconds to cut from the past end of the filter
%   endCut - seconds to cut from the future end of the filter
%   inv - 1 if response values are inverted (e.g. if you're using 
%       the voltage sensor ASAP), 0 if not
% 
% OUTPUT:
%   alignedFilters - matrix of past-future filters after each cell's filter
%       is aligned. Each column vector is the filter for one cell 
% 
% last update: 03/23/17 MX

function alignedFilters = alignFiltersXCorr(pfFilters, iRef, tau, ...
    stimIFI, startCut, endCut, inv)

% Define a time window to do the cross-correlation for any given filter
% To give enough room to shift the signals, we cut a bit of the front and
% back of the filter.
iStart = round(startCut/stimIFI); % index marking start of the short filter
iEnd = length(pfFilters(:,1)) - round(endCut/stimIFI); % index marking end
% iStart = round(length(pfFilters(:,1))/fStart); 
% iEnd = iStart*fEnd; 
alignedFilters = zeros(iEnd-iStart+1,1);

% Pick one cell as the reference signal to do cross-correlation against
ref = iRef;
refSignal = pfFilters(:,ref);
refSigWin = refSignal(iStart:iEnd); % signal within the window

% For each other cell, find the optimal delay shift for aligning that
% cell's signal against the reference signal
for r = 1:size(pfFilters,2)
    signal = pfFilters(:,r);
    sigWin = signal(iStart:iEnd);
    [corrs, lags] = xcorr(refSigWin, sigWin);
    [~,ind] = max(abs(corrs));
    delay = lags(ind); % negative = shift signal left, pos = shift right
    % shift the signal based on the delay
    alignedFilters(:,r) = signal(iStart-delay:iEnd-delay); 
end 

t = -tau:stimIFI:tau-stimIFI;
tWin = t(iStart-delay:iEnd-delay);

% plot the shifted signals on top of each other
figure;
for r = 1:length(alignedFilters(1,:))
    plot(tWin, alignedFilters(:,r)); hold on;
end 
meanPF = mean(alignedFilters,2);
plot(tWin, meanPF, 'k', 'LineWidth', 2);
line([0 0], ylim, ...
'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');
ylabel('filter amplitude(AU)'); xlabel('time prior to response');

if (inv) % reverse axes on inverted
    set(gca,'YDir','reverse');
end

end 