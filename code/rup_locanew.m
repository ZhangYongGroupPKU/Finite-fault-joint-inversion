function [loca,dep,epiloca,epidep]=rup_locanew(fault,grid,sizegrid,epi)
%==========================================================================
%[loca,dep,epiloca,epidep]=rup_loca0(fault,grid,sizegrid,epi);
% rup_loca0 is used in ruptime. it is used for calculating the positions of 
% 4 corners for a plane fault.
% 
%         strike:
%        -----------------------> 
%       /   1                 2
%      / 
%     /   3                 4           1,2,3,4 is the 4 corners 
%    dip
%
% as an example:
%   clear all
%   fault=[90,90,0];grid=[3,3];sizegrid=[100,100];epi=[36,90,2,2];
%--------------------------------------------------------------------------
%  Input:
%       fault: [stike,dip rake]
%        grid: also is the 'sizemat' in other functions.
%              num in [strike,dip,rake] for 3-d model
%              and num in [strike,dip] for 2-d model
%    sizegrid: the size in [strike,dip,rake] for 3-d model
%              and size in [strike,dip] for 2-d model. Units: km
%         epi: the first 2 elements are location of the epicenter: lat and
%              long, the next others are source position in the grids for
%              both 2-d and 3-d models
%
%  Output:
%       loca: a matrix has a size of [4,2]. first column is lat, second is
%             long. from top to bottom, are the locations of for corners:
%             1,2,3,4, respectively
%        dep: depthes for corners:1,2,3,4, respectively
%    epiloca: same as epi(1:2)
%     epidep: depth of the hypocenter
% -------------------------------------------------------------------------
%                                  Zhang Yong, Chen Yun-Tai, Xu Li-Sheng
%                                        2006/04,Peking University 
%                              repaired 2007/01/30/01:00,Peking University 
% =========================================================================

if epi(3)<0||epi(3)>grid(1)||epi(4)<0||epi(4)>grid(2)
    error('position of hypocenter should be correctly on the fault plane!')
end

strike=fault(1);
dip=fault(2).*pi./180;
fai=zeros(4,1);
f=zeros(4,1);
sita=zeros(4,1);
loca=zeros(4,2);

%calculate the angles:
if epi(3)==1
    fai(1)=pi./2;
    fai(3)=pi./2;
else
    fai(1)=atan(abs((sizegrid(2).*(epi(4)-1))./(sizegrid(1).*(epi(3)-1))));
    fai(3)=atan(abs((sizegrid(2).*(grid(2)-epi(4)))./(sizegrid(1).*(epi(3)-1))));
end
if epi(3)==grid(1)
    fai(2)=pi./2;
    fai(4)=pi./2;
else
    fai(2)=atan(abs((sizegrid(2).*(epi(4)-1))./(sizegrid(1).*(grid(1)-epi(3)))));
    fai(4)=atan(abs((sizegrid(2).*(grid(2)-epi(4)))./(sizegrid(1).*(grid(1)-epi(3)))));
end

%calculate the distances for the four conners to the epicenter
f(1)=sqrt((sizegrid(1).*(1-epi(3))).^2+(sizegrid(2).*(1-epi(4))).^2);
f(2)=sqrt((sizegrid(1).*(grid(1)-epi(3))).^2+(sizegrid(2).*(1-epi(4))).^2);
f(3)=sqrt((sizegrid(1).*(1-epi(3))).^2+(sizegrid(2).*(grid(2)-epi(4))).^2);
f(4)=sqrt((sizegrid(1).*(grid(1)-epi(3))).^2+(sizegrid(2).*(grid(2)-epi(4))).^2);

%project the distances to horizontal plane (with depth of 1 or 2)
f=f.*cos(asin(sin(dip).*sin(fai)));

%project the angles to horizontal plane (with depth of 1 or 2)
fai=atan(cos(dip).*tan(fai));

%the horizontal angles for the orientation of each distance vector
sita(1)=(strike+fai(1)*180/pi-270).*pi./180;
sita(2)=(strike-fai(2)*180/pi-270).*pi./180;
sita(3)=(strike-fai(3)*180/pi-270).*pi./180;
sita(4)=(strike+fai(4)*180/pi-270).*pi./180;

%get the positions for the 4 corner projected into the horizontal planes
loca(1,:)=f(1).*[-sin(sita(1)),cos(sita(1))./cosd(epi(1))]./111.2+epi(1:2);
loca(2,:)=f(2).*[sin(sita(2)),-cos(sita(2))./cosd(epi(1))]./111.2+epi(1:2);
loca(3,:)=f(3).*[-sin(sita(3)),cos(sita(3))./cosd(epi(1))]./111.2+epi(1:2);
loca(4,:)=f(4).*[sin(sita(4)),-cos(sita(4))./cosd(epi(1))]./111.2+epi(1:2);

%depth of 4 corners
dep1=sizegrid(2)./2;
dep2=sizegrid(2)./2;
dep3=(grid(2)-1/2).*sizegrid(2);
dep4=(grid(2)-1/2).*sizegrid(2);
dep=[dep1;dep2;dep3;dep4].*sin(dip);

%location and delpth of the hypocenter 
epiloca=epi(1:2);
epidep=(epi(4)-1/2).*sizegrid(2).*sin(dip);
return
%==================================end=====================================