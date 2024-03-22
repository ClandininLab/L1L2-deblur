%% FullFieldFlashOntoGray
% Should work whether or not you have a search/reference stimulus
% last updated: 12/5/16

%% Find the reference and fffg columns
% if roiDataMat has only one column, then that column is both the reference
% and nonref column 
if length(roiDataMat(1,:)) == 1
    refColumn = 1;
    nrefColumn = refColumn;
    singleCol = true;
else 
    singleCol = false;
    % label the columns in the data matrix
    for s = 1:length(roiDataMat(1,:))
        if strcmp(roiDataMat(1,s).stimcode, refStimCode)
            refColumn = s;
        elseif isa(roiDataMat(1,s).stimDat.obj, 'FullFieldFlashOntoGray')
            nrefColumn = s;
        end 
    end 
end 

% fetch imaging inter-frame interval
imIFI = roiDataMat(1,1).imIFI;

%% Screen responses of individual ROIs, label as responding, 
% non-responding, or inverted
if singleCol % if only one column exists
    [roiDataMat, iResp, iInv] = filterROIs_stimLockedPlots_FFFoG_noRef(roiDataMat, pairedEpochs, 1/imIFI, binWidthMult);
else % if there's more than one column
    [roiDataMat, iResp, iInv, framesPerLDCycle] = filterROIs_stimLockedPlots_FFFoG(roiDataMat, ...
    refColumn, nrefColumn, pairedEpochs, 1/imIFI, inv, yScale, binWidthMult);
end  

respROIMat = roiDataMat(iResp,:);

%% plot responses averaged over ROIs
if ~singleCol
    fffDuration = respROIMat(1,1).stimDat.obj.Duration{1}; 
    plot_FullFieldFlash(respROIMat(:,refColumn), inv, yScale, imIFI, ...
        framesPerLDCycle, fffDuration, refPlotTitle, 1);
    refFig = gcf;
end 

plot_FFFOntoGray(respROIMat(:,nrefColumn), inv, nrefPlotTitle, pairedEpochs);
nrefFig = gcf;

%% Save Data
save([savePath '/' saveName], 'roiDataMat', 'respROIMat', 'iResp', 'iInv', ...
    'binWidthMult', 'framesPerLDCycle', 'refStimCode', 'inv', 'yScale', ...
    'imIFI', '-v7.3');
% save figures
if ~singleCol
    saveas(refFig,[figPath '/' saveName '_ref'],'fig');
end 
saveas(nrefFig,[figPath '/' saveName '_fffgFig'],'fig');