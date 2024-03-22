% parseBrukerXML.m
%
% Function to parse Bruker XML file for frame times, lines per frame,
%  pixels per line. Lots of other info in the xml file that could be parsed
%  out, but ignore for now (and probably forever)
%
% Assumes consistent tree structure of Bruker XML file
%
% Note: this takes a while to run b/c parseXML() takes a while with XML
%  files of this size
%
% INPUTS:
%   xmlFilePath - full path to XML file 
%
% OUTPUTS:
%   frameTimes - vector of start times of each frame, in sec
%   xDim - number of pixels in x dimension (pixels per line)
%   yDim - number of pixels in y dimension (lines per frame)
%
% CREATED: 11/14/21 - HHY
%
% UPDATED:
%   11/16/21 - HHY
%
function [frameTimes, xDim, yDim] = parseBrukerXML(xmlFilePath)


    % read in XML file as struct
    xmlStruct = parseXML(xmlFilePath);
    
    % get names of childen in first level of XML file
    for i = 1:length(xmlStruct.Children)
        lvl1Children{i} = xmlStruct.Children(i).Name;
    end
    
    % find pixels per line info under: PVStateShard>PVStateValue
    %  key="linesPerFrame" value=""
    
    % index of PVStateShard
    pvInd = find(strcmp(lvl1Children, 'PVStateShard'));
    
    % convert PVStateShard into struct where key is field name and value is
    %  value in that field
    for i = 1:length(xmlStruct.Children(pvInd).Children)
        % not all Children have values
        if ~isempty(xmlStruct.Children(pvInd).Children(i).Attributes)
            % attributes for this Child
            thisAtt = xmlStruct.Children(pvInd).Children(i).Attributes;
            
            % if this Attribute is a struct of length 2 - others are length
            %  1 b/c they have Children; ignore these
            if (length(thisAtt) == 2)
                % get names of attributes
                for j = 1:length(thisAtt)
                    attNames{j} = thisAtt(j).Name;
                end
                % index of key
                keyInd = find(strcmp(attNames,'key'));
                % name of key
                keyName = thisAtt(keyInd).Value;
                % index of value
                valInd = find(strcmp(attNames,'value'));
                % name of value
                valName = thisAtt(valInd).Value;
                
                % add this to pvStruct
                pvStruct.(keyName) = valName;
            end
        end
    end
    
    % from pvStruct, get xDim and yDim
    yDim = str2num(pvStruct.linesPerFrame);
    xDim = str2num(pvStruct.pixelsPerLine);
    
    
    % find frame timing info under: Sequence>Frame>relativeTime (lots of
    %  nodes named "Frame"
    
    % index of Sequence
    seqInd = find(strcmp(lvl1Children, 'Sequence'));
    
    % get frame times, from all Attributes for Children of Sequence with 
    %  Name "Frame"
    
    % initialize frameTimes vector
    frameTimes = [];
    for i = 1:length(xmlStruct.Children(seqInd).Children)
        % name of this Child
        thisName = xmlStruct.Children(seqInd).Children(i).Name;
        
        % if name of this Child is "Frame", process it
        if (strcmp(thisName, 'Frame'))
            % attributes of this Child that is Frame
            thisAtt = xmlStruct.Children(seqInd).Children(i).Attributes;
            
            % get names of attributes
            for j = 1:length(thisAtt)
                attNames{j} = thisAtt(j).Name;
            end
            
            % index for relativeTime
            relTInd = find(strcmp(attNames, 'relativeTime'));
            % relative time value
            thisTime = str2num(thisAtt(relTInd).Value);
            
            % add this time to frameTimes
            frameTimes = [frameTimes; thisTime];
            
        end
    end   
end