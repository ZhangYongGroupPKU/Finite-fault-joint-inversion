function [stf,rsq,stfall,syn]=pldf(a,b,ll,n)
% =========================================================================
% stf=pldf(a,b,ll,n) --->  PLDF is to make deconlvolution by using PLD
%                          method in frequency domain
% -------------------------------------------------------------------------
% Output:
%     stf: the source time function
% Input:
%     a:   seismograms of the main shock. can be vector or a matix
%     b:   seismograms of the aftershock, as Empirical Green's Function, or
%          Calculated Green's Function. also can be a vector or a matrix
%     ll:  the length of the source time function
%     n:   the number of iteration
% -------------------------------------------------------------------------
%                     Zhang Yong
%              2006/03.   revised on  2006/10/06/18:54   Peking university
%                         revised again on  2007/06/18/21:16   µØÇòÎïÀíËù
% =========================================================================

%if nargin==3
%     n=round(pldnum_iter(b));
% end
sa=size(a);sb=size(b);
% if sb(1)>sa(1)
%     b=b(1:sa(1),:);
% end  % b can not be longer than a ,else it is not true ,so delete
len=2.^nextpow2(sa(1));
pb=fft(b,len);
t=1./max(abs(pb)).^2;
%b=[b;zeros(l-length(b),1)];

pa=fft(a,len);
xx=zeros(len,sa(2));
cpb=conj(pb);
%  prepare for the iteration
if sa(2)==1                  % while a and b are vectors
    px=t.*cpb.*pa;
    py=t.*cpb.*pb;
else                         % while a and b are matrixs
    %     px=(cpb.*pa)*diag(t);
    %     py=(cpb.*pb)*diag(t);
    px=zeros(size(pa));
    py=zeros(size(pa));
    for i=1:sa(2)
        px(:,i)=t(i).*cpb(:,i).*pa(:,i);
        py(:,i)=t(i).*cpb(:,i).*pb(:,i);
    end
end

rsq=zeros(n,sa(2));
stfall=zeros(len,sa(2),n);

%[bb,aa]=butter(3,0.2*3/2);

py=1-py;
for i=2:n                             % start iteration
    pf=fft(xx,len);
    
    z=real(ifft(pf.*pb-pa));
    rsq(i-1,:)=sum(z.*z);
    %rsq(i-1,:)=rsq(i-1,:).*rsq(i-1,:);
    
    pg=px+py.*pf;
    xx=real(ifft(pg));
    
    if isscalar(ll)
        xx(ll+1:end,:)=0;
    else
        
    end
    xx(1,:)=0;
    xx(xx<0)=0;
    
%     [xm,m]=max(xx);
%     xx(:)=0;
%     xx(m)=xm;
    
    %xx=filtfilt(bb,aa,xx);xx(xx<0)=0;
    
    
    stfall(:,:,i)=xx;
end
z=real(ifft(fft(xx,len).*pb-pa));
rsq(end,:)=sum(z.*z);

stf=zeros(len,sa(2));
for i=1:sa(2)
    [m,n]=min(rsq(:,i));
    stf(:,i)=stfall(:,i,n);
end
if nargout==4
    pstf=fft(stf);
    syn=ifft(pb.*pstf);
    syn(size(a,1)+1:end,:)=[];
end
return
%========================================================end============================================================