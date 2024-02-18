function [usedex,cita,fai]=da_del(da,inter,phase,ob)
%
if nargin==2
phase='P';
end

s=cal_slowness(da(:,1)/111.1949,0,phase);
av=111.1949./s;
if phase=='P'
    v=6;
elseif phase=='S'
    v=6/sqrt(3);
end

vav=v./av;

vav(vav>1)=1;

cita=asin(v./av)*180/pi;

fai=da(:,2);

dinter=inter(1);
ainter=inter(2);



num=0;
usedex=zeros(size(fai));
% min(cita)
% dinter
% max(cita)
if nargin<4
    for i=min(cita)-dinter/2:dinter:max(cita)+dinter/2
        ddex=cita>i&cita<i+dinter;
        for j=0:ainter:360-ainter
            adex=fai>j&fai<j+ainter;
            
            dex=ddex&adex;
            
            da0=da(dex,:);
            fdex=find(dex==1);
            if isempty(da0)
                continue;
            else
                num=num+1;
                
                mda0=sum(da0,1)./size(da0,1);
                dda=sqrt((da0(:,1)-mda0(1)).^2+(da0(:,2)-mda0(2)).^2);
                [~,dexf]=min(dda);
                usedex(num)=fdex(dexf);
            end
        end
    end
elseif nargin==4
    mob=max(abs(ob));
    for i=1:size(ob,2)
        ob(:,i)=ob(:,i)/mob(i);
    end
    for i=min(cita)-dinter/2:dinter:max(cita)+dinter/2
        ddex=cita>i&cita<i+dinter;
        for j=0:ainter:360-ainter            
            adex=fai>j&fai<j+ainter;
            dex=ddex&adex;
            
            da0=da(dex,:);
            ob0=ob(:,dex);
            fdex=find(dex==1);
            if isempty(da0)
                continue;
            else
                num=num+1;
                
                mob0=sum(ob0,2)/size(ob0,2);
                obfit=gfit0(ob0,repmat(mob0,[1,size(ob0,2)]));
                [~,dexf]=max(obfit);
                usedex(num)=fdex(dexf);
            end
        end
    end
end
usedex(num+1:end)=[];

%=================================end======================================

%dist=40;dep=0;phase='P';
function s=cal_slowness(dist,dep,phase)
%s=cal_slowness(dist,dep)
if nargin==1
    dep=0;phase='P';
end
%if phase=='P'|phase=='S'
if length(phase)==1
    t=get_time_new([0:99,(0:99)+0.01],dep,phase);
    t1=t(1:end/2);
    t2=t(end/2+1:end);
    t1(isnan(t1))=0;
    t2(isnan(t2))=0;
    s0=(t2-t1)./0.01;
    s100=interp1(0:99,s0(1:end),0:0.01:99)';
    dist_loca=round((dist)*100+1);
    s=s100(dist_loca,:);
    return
else
    %s=zeros(length(dist),length(dep));
    t1=aatime(dist,dep,phase);
    t2=aatime(dist+0.001,dep,phase);
    s=(t2-t1)./0.001;
    xs1=s<0;
    s(xs1)=0;
    xs2=s>20;
    s(xs2)=20;
    return
end
