% convertRGB2Int.m
%
% Function that takes in an [r g b] sequence in (n/3) x 3 matrix and 
%  converts it into a normalized intensity sequence of n values. This 
%  particular transform undoes convertInt2RGB() by taking the lightcrafter
%  RGB sequence of 3 6-bit patterns (G0-G5, G6-R3, R4-B1) and generating
%  the a 1D sequence of intensity values.
%
% INPUT:
%   in - RGB sequence as (n/3) x 3 matrix. RGB values between 0 and 1.
%
% OUTPUT:
%   out - 1D intensity sequence as n x 1 vector, where each [r g b] is
%       converted to 3 items in out
%
% CREATED: 8/8/19 - HHY
%
% UPDATED: 8/8/19 - HHY
%

function out = convertRGB2Int(in)
    % length of final output vector
    n = size(in, 1) * size(in, 2);
    
    % preallocate output vector
    out = zeros(n, 1);

    % undo scaling by 255, which was required for appropriate handling of
    %  RGB values by Psychtoolbox. regenerates 8-bit integer values
    rescaleIn = in * 255;
    
    % convert to binary, 8-bit
    rBin = dec2bin(rescaleIn(:,1), 8);
    gBin = dec2bin(rescaleIn(:,2), 8);
    bBin = dec2bin(rescaleIn(:,3), 8);
    
    % loop through all rgb values, do conversion back to intensity seq
    %  following 6-bit pattern of lightcrafter
    for i = 1:size(rBin, 1)
        
        % binary values for each color, flipped so least significant bits
        %  first
        gBinVal = fliplr(gBin(i,:));
        rBinVal = fliplr(rBin(i,:));
        bBinVal = fliplr(bBin(i,:));
        
        % first pattern: G0-G5
        outInd = (i-1) * 3 + 1; % index into output vector
        % extract first 6 bits of green binary, convert to integer, then
        %  convert to value b/w 0 and 1 (with 6-bit resolution) 
        out(outInd) = bin2dec(fliplr(gBinVal(1:6))) / 63;
        
        % second pattern: G6-R3
        outInd = (i-1) * 3 + 2; % index into output vector
        % extract bits 7 and 8 of green, first 4 of red, flip so least
        %  significant bit is last
        pat2Bin = fliplr([gBinVal(7:8) rBinVal(1:4)]);
        % convert to integer and then value b/w 0 and 1 (with 6-bit
        %  resolution)
        out(outInd) = bin2dec(pat2Bin) / 63;
        
        
        % third pattern: R4-B1
        outInd = (i-1) * 3 + 3; % index into output vector
        % extract bits 5-8 of green, first 2 of blue, flip so least
        %  signficant bit is last
        pat3Bin = fliplr([rBinVal(5:8) bBinVal(1:2)]);
        % convert to integer and then value b/w 0 and 1 (with 6-bit
        %  resolution)
        out(outInd) = bin2dec(pat3Bin) / 63;
    end
  
end

