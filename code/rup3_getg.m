function [g,locasub,dep]=rup3_getg(pathg,fault,grid,gridsize,source,loca,epi,srate,comdex)
%==========================================================================
% g=ids_getg(pathg,fault,grid,gridsize,source,loca,epi,srate)
% This is a function to read the Green's functions from the database for
% only use of fixing mechanism (e.g., IDS method)
%--------------------------------------------------------------------------
% Input
%    pathg: path of data holder of green's functions
%    fault: [strike, dip, rake]
%     grid: numbers of sub-faults along [dip, strike] directions
% gridsize: sub-fault size in [dip, strike] directions, unit is in km
%   source: the index of the sub-fault hypocenter located on, [dip, strike]
%     loca: location of stations,[latitude,longitude]
%      epi: epicentral location, [latitude,longitude]
%    srate: sampling rate of green's functions wanted
%
% Output
%       g: output green's functions, a 3-D matrix, the size is [leng,6,nsub,nsta]
% locasub: lcoations of sub-faults,[latitude, longitude]
%     dep: depth of sub-fault, in kilometers
%--------------------------------------------------------------------------
%             Zhang Yong, Peking University, 2014-05-06
%==========================================================================

% get the moment tensor elements from the fault geometry
%M=fp2mt_zh(1,fault(1),fault(2),fault(3));

if nargin==8
    comdex=1;
end
% transfer dip priority to strike priority
grid=grid([2,1]);gridsize=gridsize([2,1]);source=[1,source([2,1])];

% get the location and depth of sub-faults
[locasub,dep]=get_subloca(fault,grid,gridsize,source,loca,epi);

dep = dep;%ÉîÕðÐÞÕý£¡£¡£¡£¡

% get the pairs of each sub-fault and each station, and calculate the
% distances and azimuthes

nsub=size(locasub,1);
nsta=size(loca,1);
locasta=repmat_zh(loca(:,1:2),nsub);
locasuba=repmat(locasub,[nsta,1]);
depa=repmat(dep,[nsta,1]);

dasub0=da_zh(locasta,locasuba,1);
dasub=[depa,dasub0];


% read the Green's functions
g=zeros(1,6,nsub,nsta);% nsub*nsta
for i=1:6
    M=zeros(6,1);M(i)=1;
    [g0]=seekg_wang(pathg,dasub,M,srate);%depth=depth(1);
    g(1:size(g0,1),i,:)=g0(:,comdex,:);% comdex==1,2,3 for z,r,t components
end
g=reshape(g,[size(g,1),6,nsub,nsta]);
return
% %tic;g=seekg_wang(pathg,dasub,M,srate);toc;
% 
% % get the back azimuth of the subfaults relative to stations
% bda=da_zh(locasuba,locasta,1);
% % transfer the back azimuth to ray angle
% bfai=bda(:,2)+180;
% 
% %-------------------------------------------------------------------------
% % convert green's functions from (V,R,T) to (E,N,U)
% tic
% sfai=sin(bfai*pi/180);cfai=cos(bfai*pi/180);
% for i=1:size(g,3)    
% %      m=[0,0,1;sfai(i),cfai(i),0;-cfai(i),sfai(i),0];
% %      g(:,:,i)=g(:,:,i)*m;    
%      
%     v=g(:,1,i);r=g(:,2,i);t=g(:,3,i);
%     g(:,1,i)=r*sfai(i)-t*cfai(i); % E
%     g(:,2,i)=r*cfai(i)+t*sfai(i); % N
%     g(:,3,i)=v; % U
% end
% toc