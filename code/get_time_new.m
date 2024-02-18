function [t]=get_time_new(dist,dep,phase,endnar);
% ========================================================================
%function [t]=get_time(dist,dep,phase,endnar);
%
%GET_TIME is to get travel times of P or S from the given time table,which
%may be created by other applications.In this case, the time table is based
%on the IASPEI91 model.
%
% input:
%        dist: epicentral distance(s) in degree, which may be a vector with
%              the same step. Ranges 0 to 99deg
%         dep: depth of the hypocenter in kilometers, which can be a vector
%              with the same step. Ranges 0 to 99km
%       phase: phase name of interests, the option for which is 'P' or 'S'
%      endnar: a parameter which controls the input style. If ANY value is
%              given to this argument. The dist or dep above does NOT have 
%              to have the same step.  
% output:
%          t: the travel time for the given phase,epicentral distances and
%             focal depths. If nargin<4, t is a matrix with rows of length
%             (dist) and columns of length(dep); If nargin==4, the DIST and
%             the DEP should have the same size, and t will have the same
%             size as the inputs.
%--------------------------------------------------------------------------
%           Zhang Yong and Xu Lisheng, 2006/09/08, Bejing
%==========================================================================

%Load the Time_table
if nargin==2 %Default phase is P
    phase='P';
end
if upper(phase)=='P' %Load P time table
    load tp
    load tp0
    tp(tp==0)=1e9;
    tp0(tp0==0)=1e9;
    tp(1:length(tp0),:)=min(tp(1:length(tp0),:),tp0);
    tp(isnan(tp))=0;
    t0=tp;
elseif upper(phase)=='S'%Load S time table
    load ts
    load ts0
    ts(ts==0)=1e9;
    ts0(ts0==0)=1e9;
    ts(1:length(ts0),:)=min(ts(1:length(ts0),:),ts0);
    ts(isnan(ts))=0;
    t0=ts;
else
    error('The phase can not be found!');
end

if nargin<4 %Get times for the given distance and depth range

    [dist,~,xdist]=unique(dist(:));
    [dep,~,xdep]=unique(dep(:));
    
    dex_dist1=floor(dist*100)+1;
    dex_dist2=ceil(dist*100)+1;
    
    dex_dep1=floor(dep)+1;
    dex_dep2=ceil(dep)+1;
    
    wei1_dist=(ceil(dist*100)/100-dist)*100;
    wei1_dist=repmat(wei1_dist(:),[1,length(dep)]);
    
    wei2_dist=1-wei1_dist;
    
    wei1_dep=ceil(dep)-dep;
    wei1_dep=repmat(wei1_dep(:)',[length(dist),1]);
    
    wei2_dep=1-wei1_dep;
    
    t11=t0(dex_dist1,dex_dep1);
    t21=t0(dex_dist1,dex_dep2);
    t12=t0(dex_dist2,dex_dep1);
    t22=t0(dex_dist2,dex_dep2);
    
    t=t11.*wei1_dist.*wei1_dep+t12.*wei2_dist.*wei1_dep+...
        t21.*wei1_dist.*wei2_dep+t22.*wei2_dist.*wei2_dep;

    t=t(xdist,:);
    t=t(:,xdep);
else %Get times for the given pairs of distances and depths
    dex_dist1=floor(dist*100)+1;
    dex_dist2=ceil(dist*100)+1;
    
    dex_dep1=floor(dep)+1;
    dex_dep2=ceil(dep)+1;
    
    wei1_dist=(ceil(dist*100)/100-dist)*100;    
    wei2_dist=1-wei1_dist;
    
    wei1_dep=ceil(dep)-dep;
    wei2_dep=1-wei1_dep;
    
    t11=t0(dex_dist1+(dex_dep1-1)*size(t0,1));
    t21=t0(dex_dist1+(dex_dep2-1)*size(t0,1));
    t12=t0(dex_dist2+(dex_dep1-1)*size(t0,1));
    t22=t0(dex_dist2+(dex_dep2-1)*size(t0,1));
    
    t=t11.*wei1_dist.*wei1_dep+t12.*wei2_dist.*wei1_dep+...
        t21.*wei1_dist.*wei2_dep+t22.*wei2_dist.*wei2_dep;
end
return
%=====================================End==================================
