function col=jet_zh2017(num)
%=======================================================================
%  col=jet_zh(m,flag)
%-----------------------------------------------------------------------
% Input:
%       num: the number of colors in the colormap, >10 would be better
% Output:
%      col: the returned colors
%-----------------------------------------------------------------------
%                               Zhang Yong, Peking University
%=======================================================================

%num=101;
c0=jet(num);
c1=c0(1:round((num+1)/3),:);c2=c0(round((num+1)*2/3):end,:);

cadd1=interm([c0(round((num+1)/3),:);[1,1,1]],[round((num-1)/5),1]);
cadd2=interm([[1,1,1];c0(round((num+1)*2/3),:)],[round((num-1)/5),1]);

col=[c1;cadd1;cadd2;c2];

return

function J = jet(m)
%JET    Variant of HSV
%   JET(M), a variant of HSV(M), is an M-by-3 matrix containing
%   the default colormap used by CONTOUR, SURF and PCOLOR.
%   The colors begin with dark blue, range through shades of
%   blue, cyan, green, yellow and red, and end with dark red.
%   JET, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   See also HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 5.7.4.2 $  $Date: 2005/06/21 19:31:40 $

nod=4;

n = ceil(m/nod);
u = [(1:1:n)/n ones(1,n-1) (n:-1:1)/n]';
g = ceil(n/2) - (mod(m,nod)==1) + (1:length(u))';
r = g + n;
b = g - n;
g(g>m) = [];
r(r>m) = [];
b(b<1) = [];
J = zeros(m,3);
J(r,1) = u(1:length(r));
J(g,2) = u(1:length(g));
J(b,3) = u(end-length(b)+1:end);
