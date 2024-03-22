% loadROIData.m
%
% Loads pData and experimental metadata into a struct for each ROI.
%  Computes dF/F.
% Does what createROIMatrix does, but when matching ROIs across multiple
%  time series
% 
% INPUTS:
%   roiMetaMat - matching of ROIs across time series, output of a
%       matchROIsAcross[] function
%   metaDat - struct where each field is a column read out from the
%       metadata spreadsheet
%   pdatapath - path to pData folder
% 
% OUTPUTS:
%   out - Returns an R by S cell array of structs where R is the number 
%       of ROIs and S is the number of stimuli per ROI
%
% Updates
% 1/17/17 HHY - removes saving of whole stimDat in out, as that
%   makes out way too big and is unnecessary for later analysis
% 3/8/17 HHY - changes dF/F baseline definition to remove class 
%   case switching
% 3/13/17 HHY - filters 25Hz noise out of background-subtracted
%   raw fluorescence trace if dither hold wasn't on
% 3/16/17 HHY - save rawStim from obj.Out
% 3/28/17 HHY - save all fields of obj.Out except date, pdData, pdTime,
%   imFrameTime, imDataName -> more generalizable; remove saving of rawStim
%   specifically
%
function out = loadROIData(roiMetaMat, metaDat, pdatapath)


cd(pdatapath);

nROIs = size(roiMetaMat, 1);
nStim = size(roiMetaMat, 2);

% initialize RxS struct array
out = struct([]);

% for each stimulus class
for s = 1:nStim
    % for each ROI
    for r = 1:nROIs
        % fetch time series ID from metadata struct (from spreadsheet)
        out(r,s).seriesID = roiMetaMat(r,s).seriesID;
        % look up index in spreadsheet
        iMS = strcmpi(metaDat.seriesID, out(r,s).seriesID);
        
        % open the corresponding pData .mat file 
        load([out(r,s).seriesID '_pData.mat']);
        % fetch the mask index for this ROI
        imask = roiMetaMat(r,s).roiMaskInd;

        % grab relevant fields from metadata struct (from spreadsheet)
        out(r,s).flyID = metaDat.flyID(iMS);
        out(r,s).roiMask = imask;
        out(r,s).genotype = metaDat.genotype(iMS);
        out(r,s).stimcode = metaDat.stimcode{iMS};
        out(r,s).zdepth = metaDat.zdepth(iMS);
        out(r,s).ditherHold = metaDat.ditherHold(iMS);
        out(r,s).wavelength = metaDat.wavelength(iMS);
        out(r,s).pwm = metaDat.pwm(iMS); 
        
        % grab fields from the pDat struct
        out(r,s).avSignal = pDat.avSignal(imask, :);
        out(r,s).dSignal = pDat.dSignal(imask, :);
        out(r,s).nFrames = pDat.nFrames;
        out(r,s).imFrameStartTimes = pDat.imFrameStartTimes(1:length(out(r,s).dSignal));
        out(r,s).imIFI = pDat.imIFI;
        out(r,s).pStimDat = pDat.pStimDat;
%         out(r,s).stimDat = pDat.stimDat; 
        out(r,s).fileloc = pDat.fileloc;
        % save parameter values from stimDat
        for i=1:length(pDat.stimDat.obj.ParamList)
            out(r,s).stimDat.(pDat.stimDat.obj.ParamList{i}) = ...
                pDat.stimDat.obj.(pDat.stimDat.obj.ParamList{i});
        end
        % also save fields of stimDat.obj.Out, but excluding pdData,
        %  pdTime, imFrameTime, imDataName, date
        stimOutFields = fieldnames(pDat.stimDat.obj.Out);
        excludeFields = {'date','pdData', 'pdTime', 'imFrameTime', ...
            'imDataName'};
        for i = 1:length(stimOutFields)
            if ~(sum(strcmp(stimOutFields(i),excludeFields)))
                out(r,s).stimDat.(char(stimOutFields(i))) = ...
                    pDat.stimDat.obj.Out.(char(stimOutFields(i)));
            end
        end
        
        % If ROI's signal is NaN, it's a badly drawn ROI. Although this
        % should be dealt with in ROI_PD_analysis_singleTimeSeries.m,
        % for now let's set the signal to zero.
        if sum(isnan(out(r,s).dSignal)) > 0
            out(r,s).dSignal = zeros(1, length(out(r,s).dSignal));
        end 
        
        % filter out 25Hz noise if dither hold was off
        % imaging frame rate must be greater than 50Hz
        if ((~out(r,s).ditherHold) && (out(r,s).imIFI < 0.02))
            d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',24.9,'HalfPowerFrequency2',25.1, ...
               'DesignMethod','butter','SampleRate',...
               1/out(r,s).imIFI);
            out(r,s).fSignal = filtfilt(d,out(r,s).dSignal);
        else
            out(r,s).fSignal = out(r,s).dSignal;
        end
        
        % Compute dF/F
        % get baseline for F0, stimulus specific, called through stim obj
        [baselineFrameTimes, baselineSignals] = ...
            pDat.stimDat.obj.selectDFFbaseline(pDat.imFrameStartTimes,...
            out(r,s).fSignal, pDat.pStimDat.lightStartTimes, ...
            pDat.pStimDat.darkStartTimes);
        
        dFF = computeDFF(out(r,s).fSignal, out(r,s).imFrameStartTimes, ...
                   baselineSignals, baselineFrameTimes);
        out(r,s).dFF = dFF;
                      
        clear pDat 
        
    end 
end 

   
end 
