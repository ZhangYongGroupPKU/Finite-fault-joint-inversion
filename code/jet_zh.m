function col=jet_zh(m,flag)
%=======================================================================
%  col=jet_zh(m,flag)
%-----------------------------------------------------------------------
% Input:
%       m:  is a scalar to decide the number of colors in the colormap
%     flag: ==1; the deep red-red-yellow-white colormap
%           ==2; the deep red-red-yellow-green-white colormap
% Output:
%      col: the returned colors
%-----------------------------------------------------------------------
%                               Zhang Yong, CEA
%=======================================================================

num=100;
if m<num
    if flag==1
        col0=jet_zh(180,1);
    else
        col0=jet_zh(160,2);
        %col0(1,:)=[];
        col0(end,:)=[];
    end
    %col0=jet0(num,flag);
    col=zeros(m,3);
    if m<3
        error('colors should be more than 3');
    else
        col(1,:)=col0(1,:);
        col(end,:)=col0(end,:);
        for i=2:m-1
            col(i,:)=col0(round(num/(m-1)*(i-1)),:);
        end
    end
else
    col=jet0(m,flag);
end
return


function col=jet0(m,flag)
col0=jet(m);

if flag==1
    ge=round(m/8*5);
    col_2=col0(ge:m,:);
    
    ge=round(ge*3/10);
    col_1=[repmat(ones(ge-2,1),[1,2]),[1:-1/(ge-2):1/(ge-2)]'];
    
    col=[col_1;col_2];
    
else
    ge1=round(m/8*5);
    ge2=round(m/8*3);
    col=col0(ge2:m,:);
    
    num=round((ge1-ge2)/2);
    col(1:num,1)=col0(ge1:-1:ge1-num+1,1);
    
end

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
