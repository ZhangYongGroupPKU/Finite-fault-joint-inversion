%clear all
%cen=[0,0];ed=[1,1];ratio=0.5;wid=[0.02];ang=[20,50];
function [h,dot]=arrow0(cen,ed,ratio,wid,ang,linecolor,facecolor)
%==========================================================================
% [h,dot]=arrow0(cen,ed,ratio,wid,ang)
%  a function to plot a arrow
%--------------------------------------------------------------------------
%  Input
%      cen: start dot. a 2-elements vector
%       ed: end dot. a 2-elements vector
%    ratio: the ratio of arrow compared to the whole length from 'cen' to 'ed'. 0.3-0.8 is suggested.
%      wid: the width of the lines of the arrow, 0.02 is suggested.
%      ang:  a 2-elements vector to describe the arrow. the 1st element should be less than the second
%  Output
%       h: the plot handle
%     dot: the 9 dots which is used to describe the arrow 
%--------------------------------------------------------------------------
%                              Zhang Yong, 2010-02-09 10:13
%            Institute of Geophysics, China Earthquake Administration
%==========================================================================
xy=[1,1];
if nargin==5
linecolor=ones(1,3)*0.5;
facecolor=ones(1,3)*0.8;

linecolor='k';
facecolor='k';
elseif nargin==6
facecolor=ones(1,3)*0.8;
end
dist=norm(ed(:)-cen(:));
dy=ed(2)-cen(2);dx=ed(1)-cen(1);
azim=atand(dy/dx);
if dx<0;  azim=azim+180;  end


dot=zeros(9,2);

dot(1,:)=cen+dist.*wid(1).*[cosd(azim+90)*xy(1),sind(azim+90)*xy(2)];dot(8,:)=dot(1,:);
dot(2,:)=cen+dist.*wid(1).*[cosd(azim-90)*xy(1),sind(azim-90)*xy(2)];dot(9,:)=dot(2,:);

dot(5,:)=ed(:)';

cen1=cen*ratio+ed*(1-ratio);

x=dist*ratio;
l=x*tand(ang(1))*tand(ang(2))/(tand(ang(2))-tand(ang(1)));
y=l/tand(ang(2));

dot(4,:)=cen1+y.*[cosd(azim+180),sind(azim+180)]+(l-dist*wid)*[cosd(azim-90)*xy(1),sind(azim-90)*xy(2)];
dot(6,:)=cen1+y.*[cosd(azim+180),sind(azim+180)]+(l-dist*wid)*[cosd(azim+90)*xy(1),sind(azim+90)*xy(2)];

dot(3,:)=dot(2,:)+(dist*(1-ratio)-y+(l-dist*wid)/tand(ang(2))).*[cosd(azim)*xy(1),sind(azim)*xy(2)];
dot(7,:)=dot(1,:)+(dist*(1-ratio)-y+(l-dist*wid)/tand(ang(2))).*[cosd(azim)*xy(1),sind(azim)*xy(2)];

hold on

h=fill(dot(:,1),dot(:,2),facecolor,'edgecolor',linecolor,'linewidth',0.25);
%plot(dot(:,1),dot(:,2),'color',linecolor);