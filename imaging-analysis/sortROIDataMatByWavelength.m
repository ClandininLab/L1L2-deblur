% sortROIDataMatByWavelength.m
%
% Function to resort roiDataMat (output of loadROIData) such that the
%  reference wavelength is the first column for all ROIs and the other
%  wavelengths are ordered in the columns as specified. Otherwise, columns
%  don't correspond to same wavelength, but rather the order they were in 
%  the spreadsheet.
%
% INPUT:
%   roiDataMat - ROI data matrix, output of loadROIData to be sorted
%   refWavelength - wavelength of reference 
%   nRefWavelength - wavelength(s) of the non-reference time series; single
%       number or array of numbers
%
% OUTPUT:
%   refColumn - which column is reference column; currently, this always
%       outputs 1
%   sortedROIDataMat - sorted roiDataMat
%

function [refColumn, sortedROIDataMat] = sortROIDataMatByWavelength(...
    roiDataMat, refWavelength, nRefWavelength)

    % ref column is always the first column
    refColumn = 1;
    
    % copy to preserve data during sorting
    sortedROIDataMat = roiDataMat;
    
    for r = 1:size(roiDataMat,1)
        for s = 1:size(roiDataMat,2)
            currentWavelength = roiDataMat(r,s).wavelength;
            if (currentWavelength==refWavelength)
                % ref column
                sortedROIDataMat(r,refColumn) = roiDataMat(r,s);
            else
                % rest of columns are in order defined by nRefWavelength
                currentInd = find(nRefWavelength == currentWavelength);
                sortedROIDataMat(r,currentInd+1) = roiDataMat(r,s);
            end 
        end 
    end

end