% quant_shaulModeling.m
%
% quick script to generate quantification metrics (and corresponding plots)
%  of Shaul's model impulse responses
%
% Updated: 6/11/18
%

clearvars

% load data
load('/Users/hyang/Dropbox (Personal)/Clandinin Lab/Manuscripts/L1 L2 Paper/ShaulModeling/Figure4TempTraceData.mat');

numContrasts = size(L1ActCell,1);
numConditions = size(L1ActCell,2);


%% loop through all traces for L1, compute quantification metrics
for i = 1:numContrasts
    for j = 1:numConditions
        % impulse response, only 1st column
        filter = L1ActCell{i,j}(:,1)';
        % time (not equal steps for adaptation)
        t = L1TimeCell{i,j};

        % end of 2nd phase 0.25 sec from stim start
        tEndPhase2 = 0.25;
        % frame corresponding to end of second phase
        endPhase2 = find(t<tEndPhase2,1,'last');

        % framePeaks of 1st and 2nd phases
        [L1quant{i,j}.framePeak1,L1quant{i,j}.framePeak2]  = ...
            computeFramePeaks(filter,i);

        % first zero crossing
        L1quant{i,j}.frameZero1 = computeFrameZero1(filter,...
            L1quant{i,j}.framePeak1,i);

        % second zero crossing
        L1quant{i,j}.frameZero2 = computeFrameZero2(filter,...
            L1quant{i,j}.framePeak1,i);

        % tPeaks - always relative to first zero crossing
        L1quant{i,j}.tPeak1 = t(L1quant{i,j}.framePeak1 - L1quant{i,j}.frameZero1);
        L1quant{i,j}.tPeak2 = t(L1quant{i,j}.framePeak2 - L1quant{i,j}.frameZero1);

        % area under curve of first phase in units AU * sec
        %  between first and second zero crossings
        L1quant{i,j}.area1 = trapz(t(...
            L1quant{i,j}.frameZero1:L1quant{i,j}.frameZero2),filter(...
            L1quant{i,j}.frameZero1:L1quant{i,j}.frameZero2));


        % area under curve of second phase in units AU * sec
        %  between second zero crossing and endPhase2 (0.25sec from stim
        %  start)
        % dark
        L1quant{i,j}.area2 = trapz(t(...
            L1quant{i,j}.frameZero2:endPhase2),filter(...
            L1quant{i,j}.frameZero2:endPhase2));        

        % ratio of 2 area of second phase/area of first phase
        L1quant{i,j}.areaRatio = L1quant{i,j}.area2 / L1quant{i,j}.area1;

        % invert area (based on contrast) and areaRatio
        if (i == 1)
            L1quant{i,j}.area1 = L1quant{i,j}.area1 * -1;
        elseif (i == 2)
            L1quant{i,j}.area2 = L1quant{i,j}.area2 * -1;
        end
        L1quant{i,j}.areaRatio = L1quant{i,j}.areaRatio * -1;        
        
    end
end

%% loop through all traces for L2, compute quantification metrics
for i = 1:numContrasts
    for j = 1:numConditions
        % impulse response, only 1st column
        filter = L2ActCell{i,j}(:,1)';
        % time (not equal steps for adaptation)
        t = L2TimeCell{i,j};

        % end of 2nd phase 0.25 sec from stim start
        tEndPhase2 = 0.25;
        % frame corresponding to end of second phase
        endPhase2 = find(t<tEndPhase2,1,'last');

        % framePeaks of 1st and 2nd phases
        [L2quant{i,j}.framePeak1,L2quant{i,j}.framePeak2]  = ...
            computeFramePeaks(filter,i);

        % first zero crossing
        L2quant{i,j}.frameZero1 = computeFrameZero1(filter,...
            L2quant{i,j}.framePeak1,i);

        % second zero crossing
        L2quant{i,j}.frameZero2 = computeFrameZero2(filter,...
            L2quant{i,j}.framePeak1,i);

        % tPeaks - always relative to first zero crossing
        L2quant{i,j}.tPeak1 = t(L2quant{i,j}.framePeak1 - L2quant{i,j}.frameZero1);
        L2quant{i,j}.tPeak2 = t(L2quant{i,j}.framePeak2 - L2quant{i,j}.frameZero1);

        % area under curve of first phase in units AU * sec
        %  between first and second zero crossings
        L2quant{i,j}.area1 = trapz(t(...
            L2quant{i,j}.frameZero1:L2quant{i,j}.frameZero2),filter(...
            L2quant{i,j}.frameZero1:L2quant{i,j}.frameZero2));


        % area under curve of second phase in units AU * sec
        %  between second zero crossing and endPhase2 (0.25sec from stim
        %  start)
        % dark
        L2quant{i,j}.area2 = trapz(t(...
            L2quant{i,j}.frameZero2:endPhase2),filter(...
            L2quant{i,j}.frameZero2:endPhase2));        

        % ratio of 2 area of second phase/area of first phase
        L2quant{i,j}.areaRatio = L2quant{i,j}.area2 / L2quant{i,j}.area1;

        % invert area (based on contrast) and areaRatio
        if (i == 1)
            L2quant{i,j}.area1 = L2quant{i,j}.area1 * -1;
        elseif (i == 2)
            L2quant{i,j}.area2 = L2quant{i,j}.area2 * -1;
        end
        L2quant{i,j}.areaRatio = L2quant{i,j}.areaRatio * -1;        
        
    end
end

%% plot quantification (TO DO)

colormap lines
cmap = colormap;
close

fieldList = {'area1', 'area2', 'areaRatio', 'tPeak1', 'tPeak2'};
xPos = 1:8;
yScale = [0 0.05; 0 0.05; 0 2; 0 0.05; 0 0.15];

% pairings of trace quantifications to pair
whichPlot = [1 2; 1 3];

for q = 1:size(whichPlot,1)
    for m = 1:length(fieldList)
        figure;
        hold on;
        whichXPos = 1;
        % plot for L1
        for n = 1:numContrasts
            for r = 1:size(whichPlot,2)
                plot(whichXPos,L1quant{n,whichPlot(q,r)}.(fieldList{m}),...
                    'o','MarkerEdgeColor',cmap(whichXPos,:),...
                    'MarkerFaceColor',cmap(whichXPos,:),'MarkerSize',8);
                whichXPos = whichXPos + 1;
            end
        end
        % plot for L2
        for n = 1:numContrasts
            for r = 1:size(whichPlot,2)
                plot(whichXPos,L2quant{n,whichPlot(q,r)}.(fieldList{m}),...
                    'o','MarkerEdgeColor',cmap(whichXPos,:),...
                    'MarkerFaceColor',cmap(whichXPos,:),'MarkerSize',8);
                whichXPos = whichXPos + 1;
            end
        end

%         legend(fliplr(legends),'location','northeast');
        title(['Which plot ' num2str(q) ' ' fieldList{m}]);
        ylim(yScale(m,:));
        xlim([0 10]);
        set(gca,'xtick',[]);
        set(gca,'xticklabel',[]);
    end
end

%% save