% angleLawOfCosines.m
%
% Function to find angle when given 3 sides of triangle using the law of
%  cosines
%
% INPUT:
%   oppSide - length of side opposite angle to find
%   adjSide1 - length of one adjacent side
%   adjSide2 - length of the other adjacent side
%
% OUTPUT:
%   gammaAngle - in degrees, angle opposite oppSide in triangle
%

function gammaAngle = angleLawOfCosines(oppSide, adjSide1, adjSide2)
    gammaAngle = acosd((adjSide1^2+adjSide2^2-oppSide^2)/...
        (2*adjSide1*adjSide2));
end