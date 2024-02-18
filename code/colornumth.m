function [numth,inter,m]=colornumth(slip)
%==========================================================================
% [numth,inter]=colornumth(slip)
% This is a function to estimate the index for the colormap while pcolor
% the slip on the fault
%--------------------------------------------------------------------------
% Input
%         slip: slip amplitudes or the maximum slip amplitude
% Output
%        numth: number of colors of colormap
%        inter: slip amplitudes each color denoting
%            m: numth*inter-max(slip(:))
%--------------------------------------------------------------------------
%                           Zhang Yong, 2009-11-19 16:24, IGPCEA, Beijing 
%==========================================================================
%  clear all
%  slip=0.57;

mslip=max(slip(:));
numth0=repmat([6:15],[5,1]);
numth0=repmat([6:9],[5,1]);

inter0=repmat([1;2;3;4;5],[1,size(numth0,2)])*10;

rsq=zeros([size(inter0),6]);
for i=1:6
    rsq(:,:,i)=( (numth0.*inter0)-mslip*10^(i-1) ) / 10^(i-1);
end
rsq(rsq<-1e-5)=max(rsq(:));

[m,n]=min(rsq(:));

ith=ceil(n/numel(inter0));
nn=mod(n,numel(inter0));
if nn==0
    nn=numel(inter0);
end

numth=numth0(nn);inter=inter0(nn)/10^(ith-1);
