% plot_linFilters_avgImpulse.m
%
% Quick and dirty script to plot Marjorie's white noise data and impulse
% responses on top of it
%
% Updated: 9/26/17

cd('/Users/hyang/Documents/New Imaging Analysis/Marjorie/AnalyzedData/whiteNoise');
load('fff300ms_fffrand9_16ms_L2ASAP2f_refFiltered_filtered_analyzed_170915.mat');

% upSampFactor = 100;
% delay = 50; %corresponds to 6.25 ms
% shiftLinFilters = zeros(size(linFilters));
% % preshift linFilters by 6.25 ms
%     for i = 1:size(linFilters,2)
%         resp = linFilters(:,i);
%         respUpSamp = interp(resp,upSampFactor); % upsample
%         
%         % shift
%         shiftRespUpSamp = circshift(respUpSamp,-1*delay);
%         
%         % downsample
%         shiftResp = downsample(shiftRespUpSamp,upSampFactor);
%         
%         % save
%         shiftLinFilters(:,i) = shiftResp;
%     end

k = 44;

cm = colormap('lines');
close all


figure; 
hold on;
% linFiltersPlot = -1*shiftLinFilters(1:k,useROIs);
linFiltersPlot = -1*linFilters(1:k,useROIs);
for n = 1:length(useROIs)
    flyID(n) = goodROIs(n).flyID;
end

t = (ifi*-0.52):ifi:((k-1.52)*ifi);
m1 = mean(linFiltersPlot,2);
e1 = std(linFiltersPlot,[],2)./sqrt(size(linFiltersPlot,2));
h1 = plot_err_patch_v2(...
    t,...
    m1 ,...
    e1, ...
    cm(1,:),(cm(1,:)+1)/2);

% for n = 1:length(useROIs)
%     plot(0:ifi:((k-1)*ifi),linFiltersPlot(:,n));
% end
% plotting parameters
xScale = [0, (k-1)*ifi];
xlabel('time (sec)');
ylabel('normalized filter amplitude');
xlim(xScale);
ylim([-0.8 0.2]);
line([0 (k-1)*ifi],[0 0],'color',[0 0 0]); % x-axis line

title(['L2 ASAP2f white noise filter n = ' num2str(length(useROIs)) ' (' ...
    num2str(length(unique(flyID))) ')']);

% compute average impulse response
cd('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170927');
load('21Dhh_ASAP2f_200pwm.mat');

respROIMat = roiDataMat(iResp,2);

meanIR = zeros(length(respROIMat), length(respROIMat(1).rats{1}));
for n = 1:length(respROIMat)
    meanIR(n, :) = (respROIMat(n).rats{1} - respROIMat(n).rats{2})/2;
end

m2 = mean(meanIR,1);
normMean = meanIR / max(abs(m2)) * max(abs(m1));
m2 = mean(normMean,1);

e2 = std(normMean,[],1)./sqrt(size(normMean,1));
h2 = plot_err_patch_v2(...
    respROIMat(1).t{1},...
    m2 ,...
    e2, ...
    cm(2,:),(cm(2,:)+1)/2);

legend([h1 h2], 'White Noise Filter',...
    'Average Light and Dark Impulse Responses','Location','southeast');




        
