% sortROIDataMatByStim.m
%
% Function to resort roiDataMat (output of loadROIData) such that the
%  reference stimulus is the first column for all ROIs and the other
%  stimuli are ordered in the columns as specified. Otherwise, columns
%  don't correspond to same stimulus, but rather the order they were in the
%  spreadsheet.
%
% INPUT:
%   roiDataMat - ROI data matrix, output of loadROIData to be sorted
%   refStimCode - stimulus code of the reference stimulus
%   nrefStimCode - stimulus code(s) of the non-reference stimulus/stimuli;
%       should be single string or cell array of strings
%
% OUTPUT:
%   refColumn - which column is reference column; currently, this always
%       outputs 1
%   sortedROIDataMat - sorted roiDataMat
%
% UPDATED:
%   1/4/22 - HHY - rewrote to behave correctly for multiple trials of same
%   stimulus
%   3/5/22 - HHY - fixed typo in "get number of stimuli in each category 
%       of non-ref stimuli", hopefully... 

function [refColumn, sortedROIDataMat] = sortROIDataMatByStim(...
    roiDataMat,refStimCode, nrefStimCode)

    % ref column is always the first column
    refColumn = 1;
    
    % prev ind - to keep track of through loops
    prevInd = 0;
    
    % copy to preserve data during sorting
    sortedROIDataMat = roiDataMat;
    
%     for r = 1:size(roiDataMat,1)
%         for s = 1:size(roiDataMat,2)
%             currentStimulus = roiDataMat(r,s).stimcode;
%             if strcmp(currentStimulus,refStimCode)
%                 currentInd = 1; % ref column is always 1st
%                 % ref column
%                 sortedROIDataMat(r,currentInd) = roiDataMat(r,s);
%                 prevInd = currentInd;
%             % when all stimuli across trials are unique
%             elseif (size(roiDataMat,2)==(length(nrefStimCode)+1))
%                 % rest of columns are in order defined by nrefStimCode
%                 currentInd = find(strcmp(currentStimulus,nrefStimCode));
%                 sortedROIDataMat(r,currentInd+1) = roiDataMat(r,s);
%             % when there are multiple trials of the same stimulus  
%             else 
%                 thisInd = find(strcmp(currentStimulus,nrefStimCode));
%                 
%                 % if the previous stimulus was the same as this stimulus,
%                 %  update 
%                 if (thisInd == prevInd)
%                     currentInd = prevInd + 1;
%                 end
%             end 
%         end 
%     end

    % loop to get current indices for stimuli
    for r = 1:size(roiDataMat, 1)
        % initialize
        refInd = [];
        nRefInd = cell(size(nrefStimCode));
    
        % get current indices for stimuli
        for s = 1:size(roiDataMat,2)
            currentStimulus = roiDataMat(r,s).stimcode;
            
            % stimulus is reference stimulus
            if strcmp(currentStimulus, refStimCode)
                refInd = [refInd, s];
            % non-reference stimulus    
            else
                currentInd = find(strcmp(currentStimulus,nrefStimCode));
                nRefInd{currentInd} = [nRefInd{currentInd}, s];
            end
        end
        
        % get number of stimuli in each category of non-ref stimuli
        stimCnt = zeros(1,length(nRefInd)); % preallocate
        for i = 1:length(stimCnt)
            stimCnt(i) = length(nRefInd{i});
        end
        
        % convert to cumulative sum
        stimCntSum = cumsum(stimCnt);
        
        % resort in order: ref stim first, then in nrefStimCode order
        for s = 1:size(roiDataMat,2)
            if (s <= length(refInd)) % reference columns
                sortedROIDataMat(r,s) = roiDataMat(r,refInd(s));
            else % non-reference columns
                % index into nRefInd cells
                whichInd = find(s>(stimCntSum+1),1,'last');

                if isempty(whichInd)
                    whichInd = 1;
                    thisInd = s - 1;
                else
                    % index into this nRefStim
                    thisInd = s - stimCntSum(whichInd) + 1;
                end
                
                % index into roiDataMat
                thisNRefInd = nRefInd{whichInd};
                thisRoiSInd = thisNRefInd(thisInd);
                
                sortedROIDataMat(r,s) = roiDataMat(r,thisRoiSInd);
            end
        end
    end
    
end