% convertInt2RGB.m
%
% Function that takes in a normalized intensity sequence of n values to 
%  [r g b] sequence in (n/3) x 3 matrix. This particular transform is to
%  create a 1D sequence of intensity values for display on a lightcrafter 
%  at 3X the nominal frame rate by having each frame as a sequence of 3 
%  6-bit patterns (G0-G5, G6-R3, R4-B1).
%
% INPUT:
%   in - 1D normalized intensity sequence. Should be of length that is
%       divisible by 3; otherwise, extra values removed. Intensity values
%       between 0 and 1, with 0 for black and 1 for white.
%
% OUTPUT:
%   out - (n/3) x 3 [r g b] matrix, where each [r g b] has the appropriate
%       values to display 3 samples of in
%
% CREATED: 4/11/19 - HHY
%
% UPDATED:
%   4/11/19 - HHY
%   5/8/19 - HHY - wasn't using bitset correctly so green and red bits
%       weren't being set appropriately
%

function out = convertInt2RGB(in)

    lenIn = length(in); % length of input

    % number of [r g b] triples in output
    numFramesOut = floor(lenIn / 3);
    
    % preallocate
    out = zeros(numFramesOut, 3);
    
    for i = 1:numFramesOut
        r = 0;
        g = 0;
        b = 0;
        
        % Pattern 1: G0-G5 (bits 1-6 of G)
        ind = (i-1)*3 + 1; % index into in
        
        % find value of frame of in on 0 to 63 scale, by rounding
        val = round(in(ind) * 63);
        
        % get bits to flip to 1 represent this value
        gOneBits = strfind(fliplr(dec2bin(val, 6)), '1');
        
        
        % Pattern 2: G6-R3 (bits 7-8 of G, 1-4 of R)
        ind = (i-1)*3 + 2; % index into in
 
        % find value of frame of in on 0 to 63 scale, by rounding
        val = round(in(ind) * 63);
        
        % value in binary, flipped so least significant bit first
        valBin = fliplr(dec2bin(val,6)); 
        
        % get bits to flip of green; add 6 b/c it's bits 7-8 of G; append
        %   to existing bits to flip
        gOneBits = [gOneBits, (strfind(valBin(1:2), '1') + 6)];
        
        % bitset gets confused if second input is []
        if (~isempty(gOneBits))
            % flip bits of green (6 bits from pattern 1, 2 bits from 
            %   pattern 2)
            g = sum(bitset(g, gOneBits, 1, 'uint8'));
        end
        
        % get bits to flip of red (4 bits - 1 to 4), for pattern 2
        rOneBits = strfind(valBin(3:6), '1');
        
        
        % Pattern 3: R4-B1 (bits 5-8 of R, 1-2 of B)
        ind = (i-1)*3 + 3; % index into in
        
        % find value of frame of in on 0 to 63 scale, by rounding
        val = round(in(ind) * 63);
        
        % value in binary, flipped so least significant bit first
        valBin = fliplr(dec2bin(val,6)); 
        
        % get bits to flip of red (4 bits); add 4 b/c it's bits 5-8 of R;
        %   append to existing bits to flip
        rOneBits = [rOneBits, (strfind(valBin(1:4), '1') + 4)];
        
        % bitset gets confused if second input is []
        if (~isempty(rOneBits))
            % flip bits of red (4 bits from pattern 2, 4 bits from pattern
            %   3)
            r = sum(bitset(r, rOneBits, 1, 'uint8'));
        end
        
        % get bits to flip of blue (2 bits)
        bOneBits = strfind(valBin(5:6), '1');
        
        % bitset gets confused if second input is []
        if (~isempty(bOneBits))
            % flip bits of blue (2 bits - 1 and 2)
            b = sum(bitset(b, bOneBits, 1, 'uint8'));
        end
        
        % add to out; divide by 255 to convert to appropriate scale for
        %  psychtoolbox RGB display
        out(i,:) = [r g b]/255;         
    end

end