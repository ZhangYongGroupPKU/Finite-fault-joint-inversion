% clear all
% fault=[300,13,90;345,13,150;5,13,170];
% grid=[5,5;4,5;2,5];%
% %grid=[10,5;14,5;11,5];
% %grid=[15,5;14,5;12,5];
% sizegrid=[50,50;50,50;50,50];source=[1,4,2];%source=[1,9,4];%source=[1,4,2];
% epi=[3.295,95.982];
% load smdl_0530
% sta_loca=str2num(mm(:,6:25));phase='P';srate=1/5;
function [delt,loca1,dep]=ruptime_nfnew(fault,grid,sizegrid,source,sta_loca,epi,phase)
%==========================================================================
%  delt=ruptime_nf(fault,grid,sizegrid,source,sta_loca,epi,phase);
%  This is a function to calculate relative times between fault consists of
%   several segments and stations.
% -------------------------------------------------------------------------
% 
%  sgement 1
%  ----------
%            \
%             \ segment 2
%              \
%               ----------------- .......
%                  segment 3 .......
%
%   Input:
%          fault: [strike,dip]. both strike and dip are vectors with size
%                  of [num of segments,1]
%          grid:  [Nstrike, Ndip]. both strike and dip are vectors with size
%                  of [num of segments,1]
%      sizegrid:  describe the size of grid. the units are km
%        source:  a row vector, and has 3 elements. source(1) is the index
%                of the segments source locats in. source(2) is the index
%                of grid for strike direction and source(3) is the index
%                of grid for dip direction in the segment--source(1).
%      sta_loca: locations of stations. [Lat,Long], bot are vectors
%           epi: [lat,long]. the epicenter
%         phase: only can be 'P' or 'S'!
%
%   Output:
%         delt: the relative times.
%
%
%  This is a program to calculate relative times for several segments. To
%  kown more details about this, you can see ruptime which can calculate
%  relative times for 3-d fault. And for a simple plane fault, This 
%  function has the same results with ruptime.
%  
%  Two function used in this function is important: rup_locanew and get_time
% -------------------------------------------------------------------------
%                            Zhang Yong, Chen Yun-tai and Xu Li-sheng 
%                                2007/04/16/22:46  Peking University
%==========================================================================

if grid(1,2)==1
    grid(:,2)=2;
    [delt,loca1]=ruptime_nfnew(fault,grid,sizegrid,source,sta_loca,epi,phase);
    delt(2:2:end,:)=[];
    return
end

ss=size(fault);
nsub=sum(prod(grid'));

%find the location of the source, in which segment
index_s=source(1);

loca=zeros(nsub,2);
dep=zeros(nsub,1);
sumgrid=[0;cumsum(prod(grid')')];
% find location of the segment from index_s to 1:
for i=index_s:-1:1
    if i==index_s
        [loca0,dep0]=rup_locanew(fault(i,:),grid(i,:),sizegrid(i,:),[epi,source(2:3)]);
        locanext=loca0;% for next
        locasource=loca0;% for along the strike direction
        lat=reshape(loca0(:,1),[2,2])';
        long=reshape(loca0(:,2),[2,2])';
        lat=interm(lat,[grid(i,2)-1,grid(i,1)-1]);
        long=interm(long,[grid(i,2)-1,grid(i,1)-1]);
      
        loca(sumgrid(i)+1:sumgrid(i+1),1)=lat(:);
        loca(sumgrid(i)+1:sumgrid(i+1),2)=long(:);
        dep0=reshape(dep0,[2,2])';
        dep0=interm(dep0,[grid(i,2)-1,grid(i,1)-1]);
        dep(sumgrid(i)+1:sumgrid(i+1))=dep0(:);
    else
        epi_o=locanext(1,:);%[lat,long]
        [loca0,dep0]=rup_locanew(fault(i,:),grid(i,:)+[1,0],sizegrid(i,:),[epi_o,grid(i,1)+1,1]);%add a vector
        locanext=loca0;% for next
        lat=reshape(loca0(:,1),[2,2])';
        long=reshape(loca0(:,2),[2,2])';
        lat=interm(lat,[grid(i,2)-1,grid(i,1)]); %add a vector
        lat(:,end)=[]; % delete the added vector: the last one
        long=interm(long,[grid(i,2)-1,grid(i,1)]);%add a vector
        long(:,end)=[]; % delete the added vector: the last one
        loca(sumgrid(i)+1:sumgrid(i+1),1)=lat(:);
        loca(sumgrid(i)+1:sumgrid(i+1),2)=long(:);
        dep0=reshape(dep0,[2,2])';
        dep0=interm(dep0,[grid(i,2)-1,grid(i,1)-1]);
        dep(sumgrid(i)+1:sumgrid(i+1))=dep0(:);
    end
end
% find location of the segment from index_s+1 to ss(1): ss(1)---number of
% segments
if index_s<ss(1)
    for i=index_s+1:1:ss(1)
        if i==index_s+1
            epi_o=locasource(2,:);
        else
            epi_o=locanext(2,:);%[lat,long]
        end
        [loca0,dep0]=rup_locanew(fault(i,:),grid(i,:)+[1,0],sizegrid(i,:),[epi_o,1,1]);%add a vector
        locanext=loca0;% for next
        lat=reshape(loca0(:,1),[2,2])'; 
        long=reshape(loca0(:,2),[2,2])';
        lat=interm(lat,[grid(i,2)-1,grid(i,1)]); %add a vector
        lat(:,1)=[]; % delete the added vector: the first one
        long=interm(long,[grid(i,2)-1,grid(i,1)]);%add a vector
        long(:,1)=[]; % delete the added vector: the first one
        loca(sumgrid(i)+1:sumgrid(i+1),1)=lat(:);
        loca(sumgrid(i)+1:sumgrid(i+1),2)=long(:);
        dep0=reshape(dep0,[2,2])';
        dep0=interm(dep0,[grid(i,2)-1,grid(i,1)-1]);        
        dep(sumgrid(i)+1:sumgrid(i+1))=dep0(:);
    end
end

loca1=loca;
% figure
% %axis([min(loca(:,2)) max(loca(:,2)) min(loca(:,1)) max(loca(:,1))]);
% axis equal
% hold on
% sl=size(loca);
% for i=1:sl(1) 
% plot(loca(i,2),loca(i,1),'.');
% pause(0.1);
% end

% calculate the depth of the source:
dep_epi=(source(3)-1+0.5)*sizegrid(source(1),2)*sin(fault(source(1),2)*pi/180);
 

% prepare for locations between stations and epicenter:
nsta=size(sta_loca);
nsta=nsta(1);
sta_loca0=repmat(sta_loca,[nsub,1]);
loca0=repmat(loca',[nsta,1]);
loca0=reshape(loca0(:),[2,length(loca0(:))/2])';

% cal distance: sub1:nsta--->sub2:nsta--->sub3:nsta....
dist=distazim_xu(sta_loca0,loca0,1);
depth=repmat(dep',[nsta,1]);
depth=depth(:);

% use dist and depth, calculate the travel time:
delt=get_time(dist,depth,phase,1);
delt=reshape(delt,[nsta,nsub])';

% distance and travel time for epicenter:
dist_epi=distazim_xu(sta_loca,epi,1);
co_t=get_time(dist_epi,dep_epi,phase);

% relative time = delt-co_t
for i=1:nsta
    delt(:,i)=delt(:,i)-co_t(i);
end
return
%==================================end=====================================