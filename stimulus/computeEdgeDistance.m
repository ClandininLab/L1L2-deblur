% computeEdgeDistance.m
%
% function to bisect one angle of triangle to form two triangles and
% return the length of the side opposite the clip angle (some number of
% degrees of whole orignal angle)
% original triangle is defined with oppSide, adjSide, and constSide
% angle bisected is opposite constSide and the clipAngle defines the
%  triangle including adjSide and edgeDist, where clipAngle is between
%  adjSide and the unnamed side
%
% INPUTS:
%   clipAngle - in degrees, angle formed when drawing line from that vertex
%       to opposite side
%   oppSide - length of side opposite the angle being bisected
%   adjSide - length of side adjacent to clipAngle
%   constSide - length of side that is not part of new triangle defined by
%
% OUTPUTS:
%   edgeDist - length of side opposite to clipAngle in new triangle
%

function edgeDist = computeEdgeDistance(clipAngle, oppSide, adjSide, ...
    constSide)

    % angle opposite constant side; compute using Law of Cosines
    constAngle = acosd((adjSide^2+oppSide^2-constSide^2)/...
        (2*adjSide*oppSide));
    
    % last angle of newly formed triangle (other angles = clipAngle and
    %  constAngle); triangles total 180 deg
    newAngle = 180 - clipAngle - constAngle;
    
    % compute edgeDist using Law of Sines
    edgeDist = (adjSide * sind(clipAngle))/(sind(newAngle));
end