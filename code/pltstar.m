function h=pltstar(num,pos,rr,ratio,color,edcolor,wid)
%================================================================
% h=pltstar(num,pos,rr,ratio,color)
% is a function to plot stars
%----------------------------------------------------------------
% num:   type of star you select to plot, usually num>4 
% pos:   postion of the center of star
% rr:    radius of star.if numel(rr)==1, radius in x axis and y 
%        axis are of the same value, else rr(1) for radius of x axis
%        and rr(2) for that of y axis
% ratio: describe the sharpness of star
% color: facecolor of star
%edcolor:edgecolor of star
%  wid:  linewidth
%---------------------------------------------------------------
%   Zhang Yong,  2009/07/02 10:40
%===============================================================
ang=360/num*pi/180;
xx=sin(0:ang/2:2*pi-1e-4);
yy=cos(0:ang/2:2*pi-1e-4);

if numel(rr)==1
    xx=xx*rr;
    yy=yy*rr;
else
    xx=xx*rr(1);
    yy=yy*rr(2);
end
xx(2:2:end)=xx(2:2:end)*ratio;
yy(2:2:end)=yy(2:2:end)*ratio;

xx([end+1,end+2])=xx([1,2]);
yy([end+1,end+2])=yy([1,2]);

if nargin<7;    
    wid=0.5;
    if nargin<6;    
        edcolor='k';
        if nargin<5;        
            color='w';
        end
    end
end

hold on
h=fill(xx+pos(1),yy+pos(2),color,'edgecolor',edcolor,'linewidth',wid);%patch
plot(xx+pos(1),yy+pos(2),edcolor);
return