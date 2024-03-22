%% analyze_FullFieldFlash_BinaryContrast.m
% Should work whether or not you have a search/reference stimulus
% last updated: 12/5/16

%% Analyze FullFieldFlash_BinaryContrast ref stimulus responses 

% if roiDataMat has only one column, then the 'reference' column is 1
if length(roiDataMat(1,:)) == 1
    refColumn = 1;
else % find the reference column in the data matrix 
    for s = 1:length(roiDataMat(1,:))
        if strcmp(roiDataMat(1,s).stimcode, refStimCode)
            % index of the column in roiDataMat for the reference time series
            refColumn = s;
        end 
    end 
end

% fetch some parameters
imIFI = roiDataMat(1,1).imIFI; % imaging inter-frame interval
imFrameRate = floor(1/imIFI); % imaging frame rate
flashDuration = roiDataMat(1).stimDat.obj.Duration{1}; % assumes that duration of one 

% screen responses of individual ROIs, label as responding,
% non-responding, or inverted
[refSLResp, iRefResp, iRefInv, framesPerCycle, BIN_WIDTH, P_VAL_THRESH] = ...
    test_aggregate(roiDataMat(:,refColumn), flashDuration*2, 1/imIFI, ...
    inv, yScale, binWidthMult);

% filter out nonresponding cells
refROIs = refSLResp(iRefResp);
% filter out all ROIs in the data matrix based on reference ROIs 
respROIMat = roiDataMat(iRefResp,:);
% Add stimulus-locked response (SLR) and std errors to the reference column
[respROIMat(:, refColumn).SLR] = refROIs(:).rats;
[respROIMat(:, refColumn).stdErr] = refROIs(:).stdErr;

%% Plot average and also individual traces 
plot_FullFieldFlash(refROIs, inv, yScale, imIFI, framesPerCycle, flashDuration, refPlotTitle, 1);

%% Plot just the average
plot_FullFieldFlash(refROIs, inv, yScale, imIFI, framesPerCycle, flashDuration, refPlotTitle, 0);

%% Save stimulus-specific params
save([savePath '\' saveName], 'respROIMat', 'roiDataMat', 'refSLResp', ...
    'iRefResp', 'iRefInv', 'refStimCode', 'inv', ...
    'yScale', 'binWidthMult','framesPerCycle', '-v7.3');


