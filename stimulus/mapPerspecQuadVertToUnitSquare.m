% mapPerspecQuadVertToUnitSquare.m
%
% Function that takes in verticies of perspective quadrilateral (computed
%  from the function computePerspecQuadVerticies) and returns coefficients
%  of mapping
% See explanation of how to solve here: 
%  https://www.particleincell.com/2012/quad-interpolation/
%
% INPUT:
%   x - vertical vector of 4 elements, x1-x4, starting from lower left
%        corner and going around CCW; vertical dimension
%   y - vertical vector of 4 elements, y1-y4, starting from lower left
%        corner and going around CCW; horizontal dimension
%
% OUTPUT:
%   a - vertical vector of 4 elements, a1-a4, coefficients to map
%        quadrilateral to unit square, for x
%   b - vertical vector of 4 elements, b1-b4, coefficients to map
%        quadrilateral to unit square, for y
%
% Last updated: 11/27/16
%

function [a,b] = mapPerspecQuadVertToUnitSquare(x,y)
    A = [1 0 0 0; 1 1 0 0; 1 1 1 1; 1 0 1 0];
    % compute coefficients
    a = A\x;
    b = A\y;
end