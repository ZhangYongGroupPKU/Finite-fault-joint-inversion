function [t]=get_time(dist,dep,phase,endnar);
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
    tp(tp==0)=NaN;
    tp0(tp0==0)=NaN;
    tp(1:length(tp0),:)=min(tp(1:length(tp0),:),tp0);
    tp(isnan(tp))=0;
    t0=tp;
elseif upper(phase)=='S'%Load S time table
    load ts
    load ts0
    ts(ts==0)=NaN;
    ts0(ts0==0)=NaN;
    ts(1:length(ts0),:)=min(ts(1:length(ts0),:),ts0);
    ts(isnan(ts))=0;
    t0=ts;
else
    error('The phase can not be found!');
end


if nargin<4 %Get times for the given distance and depth range
    sdist=length(dist);
    sdep=length(dep);
    t=zeros(sdist,sdep);
    dist_mod=mod(dist,0.01)*100;
    dep_mod=mod(dep,1);
    for i=1:sdist
        for j=1:sdep
            if dist_mod(i)==0
                t(i,j)=t0(round(dist(i)*100)+1,fix(dep(j))+1)...
                    *(1-dep_mod(j))...
                    +t0(round(dist(i)*100)+1,fix(dep(j))+2)*(dep_mod(j));
            else
                t(i,j)=((t0(fix(dist(i)*100)+2,fix(dep(j))+1)...
                    -t0(fix(dist(i)*100)+1,fix(dep(j))+1))*dist_mod(i)...
                    +t0(fix(dist(i)*100)+1,fix(dep(j))+1))*(1-dep_mod(j))...
                    +((t0(fix(dist(i)*100)+2,fix(dep(j))+2)-t0(fix(dist(i)...
                    *100)+1,fix(dep(j))+2))*dist_mod(i)...
                    +t0(fix(dist(i)*100)+1,fix(dep(j))+2))*(dep_mod(j));
            end
            %t(i,j)=t0(round(dist(i)*100)+1,round(dep(j))+1);
        end
    end
    return
else %Get times for the given pairs of distances and depths
    sd=size(dist);
    t=zeros(sd);
    dist_mod=mod(dist,0.01)*100;
    dep_mod=mod(dep,1);
    for i=1:sd(1)
        for j=1:sd(2)
            if dist_mod(i,j)==0
                t(i,j)=t0(round(dist(i,j)*100)+1,fix(dep(i,j))+1)...
                    *(1-dep_mod(i,j))...
                    +t0(round(dist(i,j)*100)+1,fix(dep(i,j))+2)...
                    *(dep_mod(i,j));
            else
                t(i,j)=((t0(fix(dist(i,j)*100)+2,fix(dep(i,j))+1)...
                    -t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+1))*dist_mod(i,j)...
                    +t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+1))*(1-dep_mod(i,j))...
                    +((t0(fix(dist(i,j)*100)+2,fix(dep(i,j))+2)...
                    -t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+2))*dist_mod(i,j)...
                    +t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+2))*(dep_mod(i,j));
            end
            %t(i,j)=t0(round(dist(i,j)*100)+1,round(dep(i,j))+1);
        end
    end
end
return
%=====================================End==================================
