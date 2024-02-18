function [distazimp]=distazim_xu(statlocat,epilocat,flag)
%------------------------------------------------------------------------------- 
%  
%[distazimp]=distazim_xu(statlocat,epilocat); 
%    
%     DISTAZIM_XU calculates the epicentral distance and the azimuth between
%     station and source.
%Note: The input statlocat may be a matrix. In this case, output distazimp
% will be a matrix too.
%     
%     Input:
%      statlocat - [latitude longitude] of stations in degrees,'-' for 
%                   west longitude or south latitude; '+' for east
%                   longitude or north latitude
%      epilocat  - [latitude longitude] of epicenter in degrees, '-' for
%                   west longitude or south latitude; '+' for east
%                   longitude or north latitude
%      
%     Output:
%      distazimp - [delta(degree) dist(km) stataz(degree) epiaz(degree)] 
%                  of stations
%      
%     Reference: Elementary Seismology, Richter, p317ff.
%                Handbuch der Mathematik, p331
%                
%     Xu Lisheng, Beijing, Aug.2006
%------------------------------------------------------------------------------     
%Some preparations
k=pi./180;
slon=statlocat(:,2).*k;
elon=epilocat(:,2).*k;
slat=atan(0.9933.*tan(statlocat(:,1).*k));
elat=atan(0.9933.*tan(epilocat(:,1).*k));

%Calculation of distance in rad
drad=acos(sin(slat).*sin(elat)+cos(slat).*cos(elat).*cos(slon-elon));

%calculation of the azimuth
saz=acos((sin(elat)-sin(slat).*cos(drad))./(cos(slat).*sin(drad)));
eaz=acos((sin(slat)-sin(elat).*cos(drad))./(cos(elat).*sin(drad)));

% conversion rad-deg                        
delta=drad./k;
if nargin>2
    % only calculate distance 
    distazimp=delta;
    return
end

epicaz=eaz./k;
stataz=saz./k;
%distance in km
dist=delta*111.2;

% make sure the azimuth is measured from N clockwise
x=sin(slon-elon);
y=sin(elon-slon);
xx=find(x>0);
stataz(xx)=360-stataz(xx);
yy=find(y>0);
epicaz(yy)=360-epicaz(yy);
%Change variables
temp=epicaz;
epicaz=stataz;
stataz=temp;
distazimp=[delta dist stataz epicaz];
%===============================End==============================================
