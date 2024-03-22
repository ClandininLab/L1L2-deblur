% LMtoXY.m
%
% Function to convert set of (l,m) coordinates, the coordinates of the unit
%  square mapped by the perspective quadrilateral, to (x,y) coordinates,
%  the coordinates of the quadrilateral in texture degrees coordinates
% Uses definition of x and y:
%   x = a1 + a2*l + a3*m + a4*l*m
%   y = b1 + b2*l + b3*m + b4*l*m
% Reference: https://www.particleincell.com/2012/quad-interpolation/
%
% INPUTS:
%   a - vector of 4 elements, a1-a4, coefficients to map
%        quadrilateral to unit square, for x; from
%        mapPerspecQuadVertToUnitSquare function
%   b - vector of 4 elements, b1-b4, coefficients to map
%        quadrilateral to unit square, for y; from
%        mapPerspecQuadVertToUnitSquare function
%   l - input l coordinate to convert, can be single number, 1 dimensional
%        vector, or matrix; vertical dimension of unit square; must be
%        between 0 and 1, inclusive. Must be same size as m.
%   m - input m coordinate to convert, can be single number, 1 dimensional
%        vector, or matrix; horizontal dimension of unit square; must be
%        between 0 and 1, inclusive. Must be same size as l.
%
% OUTPUTS:
%   x - quadrilateral coordinates for vertical dimension; output dimensions
%        (number, vector, matrix) same as those of l
%   y - quadrilateral coordinates for horizontal dimension; output 
%        dimensions(number, vector, matrix) same as those of m
%
% Last Updated: 11/27/16
%

function [x,y] = LMtoXY(a,b,l,m)
    % error checking
    if (size(a)~=[1, 4]) & (size(a)~=[4, 1])
        disp('a must be a vector of length 4');
        x = [];
        y = [];
        return
    elseif (size(b)~=[1, 4]) & (size(b)~=[4, 1])
        disp('b must be a vector of length 4');
        x = [];
        y = [];
        return
    elseif (size(l)~=size(m))
        disp('l and m must be the same size');
        x = [];
        y = [];
        return
    end
    
    % calculate x and y
    x = a(1) + a(2)*l + a(3)*m + a(4)*(l.*m);
    y = b(1) + b(2)*l + b(3)*m + b(4)*(l.*m);
end
