function vec_out=t_shift(vec_in,sdot)
%==========================================================================
%  vec_out=t_shift(vec_in,sdot)
%  here is a function to shift the wave (vec_in) for sdot
%--------------------------------------------------------------------------
% Input
%   vec_in: the origin wave, should be a vetor or 2-dimension matrix
%   sdot: a integer scalar describing the time shift by sampling interval
%         sdot>0, pull vec_in
%         sdot<0, press vec_in
% Output
%  vec_out: shifted wave, also a vector, having the length
%           of length(vec_in)+sdot
%
% Notice: the function 't_shift' can also work as another function of
%        'resample', it is much faster than the function 'resample' of
%         Matlab, but with a slight effect of low pass filter.
%--------------------------------------------------------------------------
%        Zhang Yong, 2012-01-31 23:55, Berlin
%==========================================================================
if sdot==0
    vec_out=vec_in;
    return
end

if size(vec_in,1)==2&&sdot==-1
    vec_out=mean(vec_in);
    return;
end

dt=(size(vec_in,1)-1)/((size(vec_in,1)-1)+sdot);

dex=(1:dt:size(vec_in,1))';

% if size(vec_in,2)>1
% %dex=repmat(dex,[1,size(vec_in,2)]);
% end

onevec=ones(1,size(vec_in,2));

fdex=floor(dex);
cdex=ceil(dex);
vec1=vec_in(int32(fdex),:);
vec2=vec_in(int32(cdex),:);

%wei2=mod(dex,1);
wei2=(dex-fdex)*onevec;
wei1=1-wei2;

vec_out=vec1.*wei1+vec2.*wei2;
return
%==========================================================================
    
    
    
%The old codes: work slowly

% dt=sdot/(size(vec_in,1)-1);
% 
% dt_sh=(0:dt:sdot)';
% 
% num_dt=(1:size(vec_in,1))'+dt_sh;
% 
% num_new=(1:size(vec_in,1)+round(sdot))';
% 
% [ztemp,n]=histc(num_new,num_dt);
% 
% %temp=[num_new,[num_dt;zeros(sdot,1)],n];
% 
% vec_in(end+1)=vec_in(end);
% num_dt(end+1)=num_dt(end);
% 
% vec1=vec_in(n);
% vec2=vec_in(n+1);
% 
% wei1=num_dt(n+1)-num_new;
% wei2=num_new-num_dt(n);
% %wei=num_dt(n+1)-num_dt(n);
% wei=num_dt(2)-num_dt(1);
% wei1=wei1./wei;
% wei2=wei2./wei;
% 
% vec_out=vec1.*wei1+vec2.*wei2;
%====================================================================