% getPDDataBruker.m
%
% Function to get photodiode data from Bruker output (csv file). csv file
%  column 1 is time, in ms, and column 2 is the voltage value. Has a header
%  row.
%
% INPUTS:
%   csvFilepath - full file path to csv file of photodiode data
%
% OUTPUTS:
%   pdData - photodiode voltage values
%   pdTime - time, in sec, of each of the pdData time points
%
% CREATED: 11/14/21 - HHY
%
% UPDATED:
%   11/22/21 - HHY
%
function [pdData, pdTime] = getPDDataBruker(csvFilepath)

    % some constants
    ms2Sec = 1000; % conversion b/w ms and seconds
    numCols = 2; % number of columns in CSV file
    delim = ','; % delimiter is , in CSV
    dataFirstLine = 2; % first line of data is 2nd b/c file has header
    varNames = {'timeMs', 'pdDat'}; % variable names
    varTypes = 'double'; % data type for variables is double
    
    % create opts object for importing data
    opts = delimitedTextImportOptions('NumVariables', numCols, ...
        'VariableNames', varNames, 'Delimiter', delim, 'DataLines', ...
        dataFirstLine);
    % set variable types
    opts = setvartype(opts, varNames, varTypes);
    
    % import using readmatrix function (note: only 2019a and later)
    csvData = readmatrix(csvFilepath, opts);
    
    % pdData is column 2
    pdData = csvData(:,2);
    
    % pdTime is column 1, converted to sec from ms
    pdTime = csvData(:,1) / ms2Sec;
    
end