function [g_end,dep_out]=seekg_wang(pathdir,dda,M,rate)
%==========================================================================
% [g_end,dep_out]=seekg_wang(pathdir,dda,M,rate)
% This is a function to read the Green's functions from the database 
%--------------------------------------------------------------------------
% Input
%  pathdir: the filefolder of database
%      dda: [dep,dist,azim], azim is defined as north-east
%     depg: depth of green's function in the database
%        M: [Mxx,Mxy,Mxz,Myy,Myz,Mzz]
%     rate: the sampling rate of the output green's functions you want
%
% Output
%   g_end: output green's functions, the size of [leng,3,n_of_channels]
%--------------------------------------------------------------------------
%             Zhang Yong, Peking University, 2014-05-06
%==========================================================================
headfile=dir([pathdir,'\GreenInfo*']);

fid1=fopen([pathdir,headfile(1).name],'r');
for i=1:6;fgets(fid1);end
tpara=fscanf(fid1,'%f');

for i=1:2;fgets(fid1);end
dist0=fscanf(fid1,'%f');
dist0(1)=[];

z=fgets(fid1);
while 1
    if length(z)>25
        if strcmp(z(1:25),'#   list of source_depths')
        break;
        end
    end
    z=fgets(fid1);
end
%for i=1:20;fgets(fid1);end
numdep=str2num(fgets(fid1));
m=repmat(' ',[numdep,11]);
for i=1:numdep
    z=fgets(fid1);
    m(i,:)=z(1:11);
end
depg=str2num(m);
fclose(fid1);

leng=tpara(3);

srate=(1/tpara(2));

if nargin<4
    rate=srate;
end

%rate

dep=dda(:,1);
dist=dda(:,2);
azi=dda(:,3);

%make a transfer for azim from 'north to east' to 'south to east'
azi=180-azi;

%--------------------------------------------------------------------------
% get the timeshift of green's function
km2deg=6371*pi/180;

distdeg=dist/km2deg;dist0deg=dist0/km2deg;
%max(dep)
if max(dep)<100
    tp=get_time(distdeg,dep,'P',1); % has the same size of dist or dep
    ts=get_time(distdeg,dep,'S',1); % has the same size of dist or dep
    tp0=get_time(dist0deg,dep,'P'); % has the size [length(dist0),length(dep)]
    ts0=get_time(dist0deg,dep,'S'); % has the size [length(dist0),length(dep)]
else
    tp=dbgrn_get_time(distdeg,dep,'P',1); % has the same size of dist or dep
    ts=dbgrn_get_time(distdeg,dep,'S',1); % has the same size of dist or dep
    tp0=dbgrn_get_time(dist0deg,dep,'P'); % has the size [length(dist0),length(dep)]
    ts0=dbgrn_get_time(dist0deg,dep,'S'); % has the size [length(dist0),length(dep)]
end
%--------------------------------------------------------------------------
%max(dist(:))
%------------------------------------------------
[sdep,nd]=sort(dep);
sdist=dist(nd);
tp=tp(nd);
ts=ts(nd);
tp0=tp0(:,nd);
ts0=ts0(:,nd);

ddep=diff(sdep);
dex=[0;find(ddep>0);length(sdep)];

% for different depth
%lengg=150*srate+leng; % srate is the sampling rate of GRN FUNCs in database 
lengg=115*rate+round(leng*rate/srate); % use rate here, (rate is the sampling rate to be used)



% get the vrt components of the corresponding mechanism
% six typical source mechanisms
exp=sum(M([1,4,6]))/3;
clvd=M(6)-exp;
ss12=-M(2); %[90,90,0]; Mtp=-Mxy
ss11=(M(1)-M(4))/2; %[225,90,-180]; (Mtt-Mpp)/2=(Mxx-Myy)/2;
ds31=M(3); %[270,90,-90]; Mrt=Mxz
ds23=-M(5); %[180,90,-90]; Mpr=-Myz

sina=sin(azi*pi/180);cosa=cos(azi*pi/180);
sin2a=sin(2*azi*pi/180);cos2a=cos(2*azi*pi/180);


% new codes:
%--------------------------------------------------------------------------


g_end=zeros(lengg,3,length(dist));
dep_out=zeros(length(dist),1);
for i=1:length(dex)-1
    % dep( dex(i)+1:dex(i+1) ) has the same focal depth
    
    dep_temp=sdep(dex(i)+1);
    if sum(abs(sdep(dex(i)+1:dex(i+1))-dep_temp))>0
        error('something wrong!');
    end
    
    [~,ndep]=min(abs(depg-dep_temp));
    
    fileg=[pathdir,'\grn_d',num2str(round(depg(ndep)))];
    

    dex_temp=dex(i)+1:dex(i+1);
    
    dist_temp=sdist(dex_temp);
    tp_temp=tp(dex_temp);
    ts_temp=ts(dex_temp);
    t_temp=[tp_temp,ts_temp];
    
    t0_temp=[tp0(:,dex(i)+1),ts0(:,dex(i)+1)];% any of [dex(i)+1:dex(i+1)] is OK

    
    g_temp=getg_wang_sub0(fileg,leng,dist_temp,t_temp,t0_temp,srate,dist0,rate);
    
    %g(1:size(g_temp,1),:,dex(i)+1:dex(i+1))=g_temp; % 1
    
    gdex=nd(dex(i)+1:dex(i+1));
    lgtemp=size(g_temp,1);
    %g(1:size(g_temp,1),:,gdex)=g_temp; % equals 1+2
    

    for j=1:length(gdex)
        m1=[exp;clvd;(ss12*sin2a(gdex(j))+ss11*cos2a(gdex(j)));(ds31*cosa(gdex(j))+ds23*sina(gdex(j)))];
        m2=[(ss12*cos2a(gdex(j))-ss11*sin2a(gdex(j)));(ds31*sina(gdex(j))-ds23*cosa(gdex(j)))];
        g_end(1:lgtemp,1,gdex(j))=g_temp(:,[1,9,3,6],j)*m1;
        g_end(1:lgtemp,2,gdex(j))=g_temp(:,[2,10,4,7],j)*m1;
        g_end(1:lgtemp,3,gdex(j))=g_temp(:,[5,8],j)*m2;
        
        dep_out(gdex(j))=depg(ndep);
    end
end
%g_end=g_end(:,:,nd); % 2

return
%=======================================================================end





% old codes:
%--------------------------------------------------------------------------
g=zeros(lengg,10,length(dist));
for i=1:length(dex)-1
    % dep( dex(i)+1:dex(i+1) ) has the same focal depth
    
    dep_temp=sdep(dex(i)+1);
    if sum(abs(sdep(dex(i)+1:dex(i+1))-dep_temp))>0
        error('something wrong!');
    end
    
    [~,ndep]=min(abs(depg-dep_temp));
    fileg=[pathdir,'\grn_d',num2str(round(depg(ndep)))];
    

    dex_temp=dex(i)+1:dex(i+1);
    sdit_temp=sdist(dex_temp);
    %sazi_temp=sazi(dex_temp);
    dist_temp=sdit_temp;
    
    tp_temp=tp(dex_temp);
    ts_temp=ts(dex_temp);
    t_temp=[tp_temp,ts_temp];
    t0_temp=[tp0(:,dex(i)+1),ts0(:,dex(i)+1)];% any of [dex(i)+1:dex(i+1)] is OK
    
    g_temp=getg_wang_sub(fileg,leng,dist_temp,t_temp,t0_temp,srate,dist0,rate);
    
    %g(1:size(g_temp,1),:,dex(i)+1:dex(i+1))=g_temp; % 1
    
    g(1:size(g_temp,1),:,nd(dex(i)+1:dex(i+1)))=g_temp; % equals 1+2
end
%g(:,:,nd)=g; % 2
%------------------------------------------------

g_end=zeros(size(g,1),3,size(g,3));
for i=1:size(g,3)
%     g_end(:,1,i)=exp*g(:,1,i)+clvd*g(:,9,i)+(ss12*sin2a(i)+ss11*cos2a(i)).*g(:,3,i)+(ds31*cosa(i)+ds23*sina(i))*g(:,6,i);
%     g_end(:,2,i)=exp*g(:,2,i)+clvd*g(:,10,i)+(ss12*sin2a(i)+ss11*cos2a(i)).*g(:,4,i)+(ds31*cosa(i)+ds23*sina(i))*g(:,7,i);
%     g_end(:,3,i)=(ss12*cos2a(i)-ss11*sin2a(i))*g(:,5,i)+(ds31*sina(i)-ds23*cosa(i))*g(:,8,i);
    
    m1=[exp;clvd;(ss12*sin2a(i)+ss11*cos2a(i));(ds31*cosa(i)+ds23*sina(i))];
    m2=[(ss12*cos2a(i)-ss11*sin2a(i));(ds31*sina(i)-ds23*cosa(i))];
    g_end(:,1,i)=g(:,[1,9,3,6],i)*m1;
    g_end(:,2,i)=g(:,[2,10,4,7],i)*m1;
    g_end(:,3,i)=g(:,[5,8],i)*m2;
end
return

