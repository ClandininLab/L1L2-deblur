% collectData.m
%
% Function to collect background NIDAQ data into global variables:
%    pdData - for photodiode voltage measurements
%    pdTime - when photodiode voltage measurements were taken
%

function collectData(src, event)
    % whichScan keeps track of which scan for indexing into array
    persistent whichScan
    global pdData pdTime imFrameTime
    if isempty(whichScan)
        whichScan = 1;
    end
    lengthData = length(event.Data);
    pdData(whichScan:(whichScan+lengthData-1)) = event.Data(:,1);
    pdTime(whichScan:(whichScan+lengthData-1)) = event.TimeStamps;
    
    imFrameTime(whichScan:(whichScan+lengthData-1)) = event.Data(:,2);
    
    whichScan = whichScan + lengthData;      
end