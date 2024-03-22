% createROIMatrix.m
%
% Creates an array of structs, one for each ROI. Structs carry metadata as 
%  as dF/F computed here. For use when there is no reference stimulus.
%
%
% INPUTS:
%   tsInd - indices of particular time series from the spreadsheet 
%   metaDat - struct where each field is a column read out(r) from the
%       metadata spreadsheet
%   pdatapath - path to pData folder
% 
% OUTPUTS:
%   out(r) - a 1D array of structs where each struct contains data for a 
%       single ROI. 
%
% Updated: 2/28/17 HHY - removes saving of whole stimDat in out, as that
%  makes out way too big and is unnecessary for later analysis
% Updated: 3/8/17 HHY - incorporates MX skipping of ROIs with no signal,
%   also changes dF/F baseline definition to remove class case switching
% Updated: 3/13/17 HHY - filters 25Hz noise out of background-subtracted
%   raw fluorescence trace if dither hold wasn't on
% Updated: 3/16/17 HHY - save rawStim from obj.Out
%
function out = createROIMatrix(tsInd, metaDat, pdatapath)

cd(pdatapath);

% keep count of ROIs
r = 1;

% initialize R-length struct array 
out = struct([]);
   
for i = 1:length(tsInd)
    iMS = tsInd(i);
    % for each time series, get the series ID 
    seriesID = metaDat.seriesID{iMS};
    
    % if pData file does not exist, skip this time series
    % exists 
    if ~isfile([seriesID '_pData.mat'])
        fprintf('%s_pData.mat\n',seriesID);
        continue;
    end  
    
    % open the corresponding pData .mat file and grab fields 
    load([seriesID '_pData.mat']);
    
    % for each ROI mask in the time series, 
    for m = 1:pDat.nMasks
        
        % if this ROI doesn't have signal, skip to the next one
        if sum(isnan(pDat.dSignal(m,:))) > 0
            continue % skip to next ROI if this ROI's signal is NaN
        end 

        % grab relevant fields from master summary metadata
        out(r).seriesID = seriesID;
        out(r).flyID = metaDat.flyID(iMS);
        out(r).roiMask = m;
        out(r).genotype = metaDat.genotype(iMS);
        out(r).zdepth = metaDat.zdepth(iMS);
        out(r).ditherHold = metaDat.ditherHold(iMS);
        out(r).wavelength = metaDat.wavelength(iMS);
        
        % save info from pDat
        out(r).avSignal = pDat.avSignal(m,:);
        out(r).dSignal = pDat.dSignal(m,:);
        out(r).nFrames = pDat.nFrames;
        out(r).imFrameStartTimes = pDat.imFrameStartTimes;
        out(r).imIFI = pDat.imIFI;
        out(r).pStimDat = pDat.pStimDat;
        % save only parameter values from stimDat, as that's all that used
        %  later and it prevents stimDat.obj.Out from being saved over and
        %  over again
        for k=1:length(pDat.stimDat.obj.ParamList)
            out(r).stimDat.(pDat.stimDat.obj.ParamList{k}) = ...
                pDat.stimDat.obj.(pDat.stimDat.obj.ParamList{k});
        end 
        % also save rawStim from obj
        if (isfield(pDat.stimDat.obj.Out,'rawStim'))
            out(r).stimDat.rawStim = pDat.stimDat.obj.Out.rawStim;
        end
        % also save epochs presented (for naturalistic stimulus)
        if (isfield(pDat.stimDat.obj.Out,'epochsPresented'))
            out(r).stimDat.epochsPresented = ...
                pDat.stimDat.obj.Out.epochsPresented;
        end
        
        % filter out 25Hz noise if dither hold was off
        % imaging frame rate must be greater than 50Hz
        if ((~out(r).ditherHold) && (1/out(r).imIFI > 50))
            d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',24.9,'HalfPowerFrequency2',25.1, ...
               'DesignMethod','butter','SampleRate',...
               1/out(r).imIFI);
            out(r).fSignal = filtfilt(d,out(r).dSignal);
        else
            out(r).fSignal = out(r).dSignal;
        end
        

        % Compute dF/F
        % get baseline for F0, stimulus specific
        [selectFrameTimes, selectSignals] = ...
            pDat.stimDat.obj.selectDFFbaseline(pDat.imFrameStartTimes,...
            out(r).fSignal, pDat.pStimDat.lightStartTimes, ...
            pDat.pStimDat.darkStartTimes);
        
        % compute dF/F
        dFF = computeDFF(out(r).fSignal, out(r).imFrameStartTimes, ...
            selectSignals, selectFrameTimes);
        out(r).dFF = dFF;
        
        % update the roi counter 
        r = r + 1;        
    end 
    
end 

out = out'; % column vector

clear pDat

end 