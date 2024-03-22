% computePerspecQuadVerticies.m
%
% Function that takes in measurements of the fly to the screen as well as
%  the angular size of a square texture that the perspective quadrilateral
%  will fall within and returns the verticies of the quadrilateral as
%  coordinates with units degrees. Texture spans units (0,0) to
%  (size,size). x is vertical dimension, y is horizontal dimension. Assumes
%  that the upper right vertex of quadrilateral matches upper right corner
%  of texture (i.e. upper right corner of screen is closest to fly). Fixes
%  upper left vertex of quadrilateral using difference in extent spanned by
%  left and right sides from fly horizontal to horizontal of top of screen.
%
% INPUTS:
%   ul - fly to upper left corner of screen
%   ur - fly to upper right corner of screen
%   ll - fly to lower left corner of screen
%   lr - fly to lower right corner of screen
%   x - vertical extent of screen
%   y - horizontal extent of screen
%   hl - fly horizontally out to vertical plane of screen, left side
%   hr - fly horizontally out to vertical plane of screen, right side
%   vl - vertical distance between horizontal plane of fly and screen, left
%         side
%   vr - vertical distance between horizontal plane of fly and screen,
%         right
%   txtrX - x coordinate, in degrees, of upper right corner of texture
%   txtrY - y coordinate, in degrees, of upper right corner of texture
% 
% OUTPUTS:
%   urx - upper right vertex, x coordinate
%   ury - upper right vertex, y coordinate
%   ulx - upper left vertex, x coordinate
%   uly - upper left vertex, y coordinate
%   llx - lower left vertex, x coordinate
%   lly - lower left vertex, y coordinate
%   lrx - lower right vertex, x coordinate
%   lry - lower right vertex, y coordinate
%
% Last Updated: 1/15/17
%

function [urx, ury, ulx, uly, llx, lly, lrx, lry] = ...
    computePerspecQuadVerticies(ul, ur, ll, lr, x, y, hl, hr, vl, vr, ...
    txtrX, txtrY)
    % error checking - return if upper right corner isn't closest to fly
    if ((ur > ul)||(hr > hl))
        disp('The upper right corner must be the one closest to the fly');
        urx = [];
        ury = [];
        ulx = [];
        uly = [];
        llx = [];
        lly = [];
        lrx = [];
        lry = [];
        return
    end
    
    % output quadrilateral upper right vertex (same as upper right corner
    %  of texture)
    urx = txtrX;
    ury = txtrY;
    
    % use Law of Cosines to compute sides of quadrilateral, in degrees
    ul_y_ur = angleLawOfCosines(y,ul,ur);
    ul_x_ll = angleLawOfCosines(x,ul,ll);
    ll_y_lr = angleLawOfCosines(y,ll,lr);
    ur_x_lr = angleLawOfCosines(x,ur,lr);
    
    % use Law of Cosines to compute diagonals of quadrilateral, in degrees
    % first compute diagonal of screen, in cm
    d = sqrt(x^2+y^2);
    % diagonals
    ul_d_lr = angleLawOfCosines(d,ul,lr);
    ur_d_ll = angleLawOfCosines(d,ur,ll);
    
    % use Law of Cosines to compute vertical angle spanned above screen on
    %  left and right sides, in degrees
    hr_vr_ur = angleLawOfCosines(vr,hr,ur);
    hl_vl_ul = angleLawOfCosines(vl,hl,ul);    
    
    % upper left vertex, x coordinate = txtrX - (hr_vr_ur-hl_vl_ul)
    ulx = txtrX - (hr_vr_ur-hl_vl_ul);
    % upper left vertex, y coordinate = txtrY - use pythagorean theorem to
    %  compute from ul_y_ur and (hr_vr_ur-hl_vl_ul); then subtract from
    %  txtrY
    uly = txtrY - sqrt(ul_y_ur^2-(hr_vr_ur-hl_vl_ul)^2);
    
    % lower left vertex
    %  compute by solving for 3 angles that divide vertical line that goes
    %  through upper left vertex - 2 angles can be computed from triangles
    %  using the law of cosines
    % upper triangle
    gamma1 = angleLawOfCosines(txtrY-uly, hr_vr_ur-hl_vl_ul, ul_y_ur);
    % quadrilateral triangle determined by one diagonal
    gamma2 = angleLawOfCosines(ur_d_ll, ul_x_ll, ul_y_ur);
    % line is 180 degrees
    gamma3 = 180 - gamma1 - gamma2;
    % gamma3 as well as ul_x_ll are angle and side of right triangle formed
    %  between upper left vertex and lower left vertex; use definition of
    %  sine and cosine to get orthogonal sides, used for (x,y)
    lly = uly + (ul_x_ll * sind(gamma3));
    llx = ulx - (ul_x_ll * cosd(gamma3));
    
    % lower right vertex
    %  compute by solving for 3 angles that divide upper right corner of
    %  texture - 2 angles determined by triangles
    % upper triangle
    gamma4 = 90 - gamma1; % same triangle as for lower left vertex
    % quadrilateral triangle determined by other diagonal
    gamma5 = angleLawOfCosines(ul_d_lr, ul_y_ur, ur_x_lr);
    % upper right corner is 90 degrees
    gamma6 = 90 - gamma4 - gamma5;
    % gamma6 as well as ur_x_lr are angle and side of right triangle formed
    % between lower right vertex and upper right vertex; use definition of
    % sine and cosine to get orthogonal sides, used for (x,y)
    lry = txtrY - (ur_x_lr * sind(gamma6));
    lrx = txtrX - (ur_x_lr * cosd(gamma6));
end

