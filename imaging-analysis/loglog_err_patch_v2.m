% loglog_err_patch_v2.m
%
% Function to plot trace (as line) with shading behind it for error values
% On semilog scale for y axis
%
% INPUT:
%   x - x values
%   y - y values for line
%   e - error values, plotted with shading as +/- that error value in y
%       dimension at corresponding x value
%   col1 - color of line
%   col2 - color of shading
%   style - linestyle parameter of Matlab plot function, optional argument;
%       changes line
%
% OUTPUT:
%   g - handle to line
%
function g=loglog_err_patch_v2(x,y,e,col1,col2,style)

x=x(:)';
y=y(:)';
e=e(:)';

ye1=y+e;
ye2=y-e; ye2=ye2(end:-1:1);
ye=[ye1,ye2];
xe=[x,x(end:-1:1)];

if(nargin==6)
    g=loglog(x,y,'color',col1,'linewidth',2,'linestyle',style);
else
    g=loglog(x,y,'color',col1,'linewidth',2);
end

hold on;
h=patch(xe,ye,col2,'linestyle','none');
uistack(g);
