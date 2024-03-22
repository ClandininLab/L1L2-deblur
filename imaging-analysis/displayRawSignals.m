% displayRawSignals.m
%
% Function that plots dF/F over time for each ROI, as well as raw
%  fluorescence before background subtraction as separate plot. Displays
%  average image with ROIs overlaid
%
% INPUT:
%   in - pData struct created by going through
%       ROI_PD_analysis_singleTimeSeries.m
%
% OUTPUT:
%   none but generates figures as side effects
%
% Last updated: 3/7/17 HHY (removed output, as it wasn't doing anything)
%

function displayRawSignals(in)

nframes = in.nFrames;
fps = 1/in.imIFI; % frames per second
nMasks = in.nMasks;

% % logical for when stimulus changes, -1, 0, +1 - for direction of change
% thresh = 0.5;
% mask = out.fstimval>thresh;
% dmask = mask(2:end)-mask(1:end-1);
% dmask = [0; dmask];


% plot signals

cm = colormap('lines');


% plots intensity over time, no background subtraction
% figure;
for i = 1:nMasks
    if (i>64)
        cmInd = mod(i,65)+1;
    else
        cmInd = i;
    end
    plot((1:nframes)/fps,in.avSignal(i,:),'color',cm(cmInd,:)); hold on;
end
xlabel('time (sec)');
title('Intensity in ROIs - before background substraction');

% plot dF/F
figure; 
for i = 1:nMasks
    if (i>64)
        cmInd = mod(i,65)+1;
    else
        cmInd = i;
    end
    % ***NOTE: this dF/F calculation is for visualization purposes. It
    % assumes the baseline F is the mean intensity of the bkgnd-subtracted 
    % average image. Depending on the stimulus, the baseline F may not be 
    % the mean! 
    
    %MMP 230301 changed while debugging fake natstim data that had NaNs
%    trace = in.dSignal(i,:)/mean(in.dSignal(i,:)) + 1*(i-1); % dF/F for each ROI
    trace = in.dSignal(i,:)/mean(in.dSignal(i,:),'omitnan') + 1*(i-1); % dF/F for each ROI
    plot((1:nframes)/fps, trace, 'color', cm(cmInd,:), 'linewidth', 2);
    hold on;
end

% % diagrams stimulus using stimulus value at each imaging frame
plot((1:nframes)/fps, in.pStimDat.stimvalIF(1:nframes)*0.2,'LineWidth',2);

axis([0 nframes/fps 0 i+1]);

% % dotted line for every stimulus transition
% inds = find(dmask~=0);
% for k = 1:length(inds)
%     if(dmask(inds(k))>0)
%         line([inds(k)/fps inds(k)/fps],[0 nMasks+1],'k','LineWidth',1,'LineStyle','-');
%     else
%         line([inds(k)/fps inds(k)/fps],[0 nMasks+1],'k','LineWidth',1,'LineStyle','--');
%     end
% end
xlabel('time (sec)'); 
ylabel('ROI #');
title('dF/F');

% % if number of imaging frames doesn't equal number of imaging frames as
% %  defined by stimulus file
% if(nframes~=length(out.fstimval))
%   warning=['Warning the response is off from the correct number of frames by: ' num2str(length(out.fstimval)-nframes)];
%   disp(warning); 
% end

% Create a colored map of ROIs
CMask = zeros(in.yPixels ,in.xPixels, 3);
for i = 1:nMasks
    if (i>64)
        cmInd = mod(i,65)+1;
    else
        cmInd = i;
    end
    curColor = cm(cmInd,:);
    curMask = cat(3,curColor(1).*in.roiMasks{i},curColor(2).*in.roiMasks{i},curColor(3).*in.roiMasks{i});
    CMask = CMask + curMask;
end
% if(~isfield(in,'AV'))
%     AV = squeeze(sum(in.alSeries,3))/nframes; % The average image
%     AV = im2double(AV);
%     AV = AV./max(AV(:));
% end
figure;imshow(in.avgIm,[]); 
hold on;h = imshow(CMask);
set(h,'AlphaData',0.5);
title('ROI masks');

end