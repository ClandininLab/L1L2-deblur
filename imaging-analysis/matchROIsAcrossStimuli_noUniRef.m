% matchROIsAcrossStimuli_noUniRef.m
%
% Function to match ROIs from the same imaging region taken over multiple
%   time series
% For when there's no unique reference stimulus
% Approximate equivalent to create_neuron_structure_all, only it matches
%   ROIs across time series
%
% INPUT:
%   inds - indicies into spreadsheet of time series to consider
%   refStimCode - stimulus code of reference stimulus (match ROIs to ROIs
%       in this time series)
%   metaDat - metadata from spreadsheet, as struct, from
%       readDatabase_selectSamples
%
% OUTPUT:
%   out - R by S cell matrix where R is the total number of ROIs and S is
%       the number of time series for a single ROI
%

function out = matchROIsAcrossStimuli_noUniRef(tsInds, refStimCode, metaDat)

    % of the set of selected time series, find the indices of the series
    % that presented the reference stimulus
    % which pData references of inds presented the reference stimulus
    refStimInd = find(strcmpi(refStimCode, metaDat.stimcode(tsInds)));

    % parameters of the selected time series 
    % limit these to the ones referenced by inds
    seriesID = metaDat.seriesID(tsInds);
    zdepth = metaDat.zdepth(tsInds);
    flyID = metaDat.flyID(tsInds);
    
    % counter for all cells in all reference stimuli pData files
    roicount = 1;
    % initialize output cell array
    out = struct([]);

    % previous zdepth and fly id, only use first of stim as reference
    prevZ = -1;
    prevFlyID = 0;

    % loop over all reference stimuli instances, match ROIs to each of
    % these
    for iRef = 1:length(refStimInd)
        thisZ = zdepth(refStimInd(iRef));
        thisFlyID = flyID(refStimInd(iRef));
        thisSeriesID = seriesID(refStimInd(iRef));

        % has to be different zdepth (or fly), don't repeat
        if ~((prevZ == thisZ) && (prevFlyID == thisFlyID))
    
            % pData file names of all time series that imaged the same cells
            sameRegionIDs = seriesID(logical((zdepth == thisZ).*(flyID == thisFlyID)));
            % pData file names without reference stimulus pData file
            sameRegionOtherStimIDs = sameRegionIDs(~strcmpi(thisSeriesID,sameRegionIDs));
            
            % if there is no nonref stimulus associated with this reference
            % time series, go to the next reference times series
            if isempty(sameRegionOtherStimIDs)
                continue
            end 
    
            % if pData for reference stim doesn't exist, skip this time series
            if ~exist([char(thisSeriesID) '_pData.mat'], 'file')
                continue;
            end 
            
            % load pData for reference stim
            load([char(thisSeriesID) '_pData.mat']);
            refStimPData = pDat;
            % load for non-reference stim
            nonrefPData = [];
            % debug
    %         sprintf('length sameRegionOtherStim = %d', length(sameRegionOtherStim))
            for iNonRef = 1:length(sameRegionOtherStimIDs)
                % if pdata file doesn't exist, skip to next iteration
                if ~exist([char(sameRegionOtherStimIDs(iNonRef)) '_pData.mat'], 'file')
                    continue
                end 
                
                % load pdata file 
                load([char(sameRegionOtherStimIDs(iNonRef)) '_pData.mat']);
                thisPData = pDat;
                nonrefPData = [nonrefPData; thisPData];
            end
    
            numRows = length(nonrefPData)+1;
            % display average images with and without masks to allow comparison
            cm = colormap('lines');
            close
    
            figure
    
            % get masks for reference stimulus
            CMask = zeros(refStimPData.yPixels, refStimPData.xPixels, 3);
            for k = 1:refStimPData.nMasks
                curColor = cm(k,:);
                curMask = cat(3,curColor(1).*refStimPData.roiMasks{k},...
                    curColor(2).*refStimPData.roiMasks{k},...
                    curColor(3).*refStimPData.roiMasks{k});
                CMask = CMask + curMask;
            end
            
            % display average image with masks
            subplot(numRows,2,1)
            imshow(refStimPData.avgIm,[]);
            hold on;
            h = imshow(CMask);
            set(h,'AlphaData',0.5); % make color mask slightly transparent
            
            % display average image, no masks
            subplot(numRows,2,2)
            imshow(refStimPData.avgIm,[]);
            
            % display average image, with and without masks, for other stimuli
            for iNonRef = 1:length(nonrefPData)
                CMask = zeros(nonrefPData(iNonRef).yPixels,...
                    nonrefPData(iNonRef).xPixels, 3);
                
                for m = 1:nonrefPData(iNonRef).nMasks
                    curColor = cm(m,:);
                    curMask = cat(3,curColor(1).*nonrefPData(iNonRef).roiMasks{m},...
                        curColor(2).*nonrefPData(iNonRef).roiMasks{m},...
                        curColor(3).*nonrefPData(iNonRef).roiMasks{m});
                    CMask = CMask + curMask;
                end
                subplot(numRows, 2, iNonRef*2+1)
                imshow(nonrefPData(iNonRef).avgIm,[]);
                hold on;
                h = imshow(CMask);
                set(h,'AlphaData',0.5);
                subplot(numRows, 2, iNonRef*2+2)
                imshow(nonrefPData(iNonRef).avgIm,[]);
            end
    
            % get matching ROIs
            for n = 1:refStimPData.nMasks
                prompt = sprintf('Enter the corresponding indicies for ROI %d: ', n);
                userInput = input(prompt);
                
                % debug
    %             sprintf('length userInput = %d', length(userInput))
    %             sprintf('length nonrefPdata = %d', length(nonrefPData))
                % input should have as many entries as time series
                if(length(userInput) == length(nonrefPData))
                    % no matching if ROI can't be found in all time series (no
                    %  match indicated with -1)
                    if(userInput > 0)
    %                     refStrct.seriesID = char(thisSeriesID);
    %                     refStrct.roiMaskInd = n;
    %                     refStrct.stimclass = class(refStimPData.stimDat.obj);
    %                     out(roicount, 1) = refStrct;
    
                        out(roicount, 1).seriesID = char(thisSeriesID);
                        out(roicount, 1).roiMaskInd = n;
                        out(roicount, 1).stimclass = class(refStimPData.stimDat.obj);
                        
                        for o = 1:length(userInput)
    %                         nonRefStrct.seriesID = char(sameRegionOtherStimIDs(o));
    %                         nonRefStrct.roiMaskInd = userInput(o);
    %                         nonRefStrct.stimclass = class(nonrefPData.stimDat.obj);
    %                         out(roicount, o+1) = nonRefStrct;
    
                            out(roicount, o+1).seriesID = char(sameRegionOtherStimIDs(o));
                            out(roicount, o+1).roiMaskInd = userInput(o);
                            out(roicount, o+1).stimclass = ...
                                class(nonrefPData(o).stimDat.obj);
                        end
    
                        roicount = roicount + 1;
                    end
                end
            end
            prevZ = thisZ;
            prevFlyID = thisFlyID;
        end
    end
    close all
end

% DO THIS CONVERSION: cell2mat

