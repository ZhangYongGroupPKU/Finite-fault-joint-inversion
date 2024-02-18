function [cen,ruptime]=astf_analysis(rtime,latlim,lonlim,num,deplim,phase,epi,dep0,loca)
%==========================================================================
%[cen,ruptime]=astf_analysis(rtime,latlim,lonlim,deplim,num,phase,epi,dep0,loca)
%--------------------------------------------------------------------------
% Input
%   rtime: the time of sub-events, in second
%  latlim: latitude range
%  lonlim: longitude range
%  deplim: depth range
%     num: number of potential location in horizontal
%  deplim: depth range, depth will vary as deplim(1):deplim(2)
%   phase: 'P' or 'S'
%     epi: epicentre
%    dep0: hypocentral depth
%    loca: locations of sations
%
% Output
%    cen: location of sub-event, [lat,long,dep]
%ruptime: rupture time of the sub-event, in second
%--------------------------------------------------------------------------
%              Zhang Yong, Peking University, 2014-07-03 22:33
%==========================================================================
% clear all
% close all
% load astf
% load data1
% [m,n]=max(stf);
% azi=da(:,2);
% rtime=n-1;

% num=50;
% lonlim=[-71,-69];
% latlim=[-20.5,-18.5];

da=da_zh(loca,epi,2);

[xx,yy,zz]=meshgrid(lonlim(1):diff(lonlim)/(num-1):lonlim(2),latlim(1):diff(latlim)...
    /(num-1):latlim(2),deplim(1):deplim(2));

locasub=[yy(:),xx(:)];
locasub=repmat(locasub,[nsta,1]);
zz=repmat(zz(:),[nsta,1]);
loca=repmat_zh(loca,numel(xx));

dasub=da_zh(locasub,loca,2);
t=get_time_new(dasub(:,1)/111.1949,zz(:),phase,1);
t=reshape(t,[numel(xx),nsta]);

t0=get_time_new(da(:,1)/111.1949,dep0,phase);
t0=repmat(t0,[1,numel(xx)]);
t0=t0';

t1=t-t0;

rtime=repmat(rtime(:)',[numel(xx),1]);


co=zeros(mean(rtime(1,:))*2,2);
for tr=1:mean(rtime(1,:))*2
    t2=tr+t1;
    [m,n]=min(sum(abs(rtime'-t2').^1));
    co(tr,:)=[m,n];
end
[~,nn]=min(co(:,1));
cen=[locasub(co(nn,2),:),zz(co(nn,2))];
ruptime=nn;

% figure
% plot(azi,[t1(co(nn,2),:)'+nn,rtime(1,:)'],'o')
% 
% 
% figure
% plot(epi(2),epi(1),'kp');hold on
% plot(cen(2),cen(1),'ks')
