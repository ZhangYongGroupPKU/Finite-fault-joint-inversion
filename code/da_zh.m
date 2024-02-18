function da=da_zh(loca,epi,flag)
%==========================================================================
%  da=da_zh(loca,epi,flag)
%  to calculate the distance and azimuth
%--------------------------------------------------------------------------
%  Input
%    loca: [lat,long], the locations of stations
%     epi: epicenter
%   if nargin<3, calculate with ellipse model, else with spherical model
% Output
%      da: distance and azimuth,[distance, azimuth]
%--------------------------------------------------------------------------
%       Zhang Yong, 2012-02-09 12:49, GFZ, Potsdam
%==========================================================================
if nargin<3
    % Supply WGS84 earth ellipsoid axis lengths in kilometers:
    a = 6378.137; % definitionally
    b = 6356.75231424518; % computed from WGS84 earth flattening coefficient definition
    c = sqrt(a*a-b*b);
    e=c./a;
    
    [dist,azi]=distance('gc',epi,loca,[a,e]);
    %[d,a]=distance(epi,loca,[a,e]);
    da=[dist,azi];
elseif flag==1
    [dist,azi]=distance(epi,loca);
    da=[dist*6371*pi/180,azi];
else
    k=pi./180;
    slon=loca(:,2).*k;
    elon=epi(:,2).*k;
    % slat=atan(0.9933.*tan(statlocat(:,1).*k));
    % elat=atan(0.9933.*tan(epilocat(:,1).*k));
    slat=atan(tan(loca(:,1).*k));
    elat=atan(tan(epi(:,1).*k));
    
    %Calculation of distance in rad
    drad=acos(sin(slat).*sin(elat)+cos(slat).*cos(elat).*cos(slon-elon));
    
    %calculation of the azimuth
    saz=acos((sin(elat)-sin(slat).*cos(drad))./(cos(slat).*sin(drad)));
    eaz=acos((sin(slat)-sin(elat).*cos(drad))./(cos(elat).*sin(drad)));
    
    % conversion rad-deg
    delta=drad./k;
    
    epicaz=real(eaz./k);
    stataz=real(saz./k);
    %distance in km
    dist=delta*6371*pi/180;
    
    % make sure the azimuth is measured from N clockwise
    x=sin(slon-elon);
    y=sin(elon-slon);
    xx=(x>0);
    stataz(xx)=360-stataz(xx);
    yy=(y>0);
    epicaz(yy)=360-epicaz(yy);
    %Change variables
    temp=epicaz;
    epicaz=stataz;
    stataz=temp;
    da=[dist,stataz,epicaz];
end
return

