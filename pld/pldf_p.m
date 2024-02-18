function [stf,rsq,stfall]=pldf_p(obs,g,ll,itern,aerfa)
% =========================================================================
% stf=pldf(a,b,ll,n) --->  PLDF_P is to make deconlvolution by using PLD
%                          method with preconditioning in frequency domain
% -------------------------------------------------------------------------
% Output:
%     stf: the source time function
% Input:
%     a:   seismograms of the main shock. can be vector or a matix 
%     b:   seismograms of the aftershock, as Empirical Green's Function, or
%          Calculated Green's Function. also can be a vector or a matrix
%     ll:  the length of the source time function
%  itern:  the number of iteration
%  aerfa:  the damp parameter. it is a small number: such as 0.01,0.02,...
%  Please see the paper of M.Bertero et al.(1997)
% -------------------------------------------------------------------------
%                    Zhang Yong , Xu Lisheng and Chen Yuntai
%                                       2007/06/18/22:11   地球物理所
% =========================================================================
so=size(obs);
sg=size(g);

if so(1)<sg(1)
    error('length of green function should be shorter than data!!');
end
if so(2)~=sg(2)
    error('number of !!');
end


fftlen=2^nextpow2(so(1));
pobs=fft(obs,fftlen);
pg=fft(g,fftlen);
pgg=conj(pg);

ll=min(ll,fftlen);

Mpg=pgg.*pg; %size: fftlen*so(2)
Mapg=max(Mpg);%size: 1*so(2)
%D=1./(Mpg+aerfa.*Mapg);%size: fftlen*so(2)

to=zeros(so(2),1);
M1=zeros(fftlen,so(2));
M2=zeros(fftlen,so(2));
D=zeros(fftlen,so(2));
for i=1:so(2)
    D(:,i)=1./(Mpg(:,i)+aerfa.*Mapg(i));
    to(i)=1+aerfa;%min((Mpg(:,i)+aerfa.*Mapg(i))./Mpg(:,i));
    M1(:,i)=to(i).*D(:,i).*Mpg(:,i);%M1=1-to.*D.*G'*G
    M2(:,i)=to(i).*D(:,i).*pgg(:,i).*pobs(:,i);%M2=to.*D.*G'*obs
end
M1=1-M1;

stfall=zeros(fftlen,so(2),itern);
rsq=zeros(itern,so(2));
xx=zeros(fftlen,so(2));
for i=2:itern                             % start iteration
    pf=fft(xx,fftlen);
    
    z=real(ifft(pf.*pg-pobs));
    rsq(i-1,:)=sum(z.*z);
    %rsq(i-1,:)=rsq(i-1,:).*rsq(i-1,:);
        
    pf1=M2+M1.*pf;
    xx=real(ifft(pf1));
    
    xx(ll+1:end,:)=0;
    xx(1,:)=0;
    xx(xx<0)=0;
    
    stfall(:,:,i)=xx;
end
z=real(ifft(fft(xx,fftlen).*pg-pobs));
rsq(end,:)=sum(z.*z);

stf=zeros(fftlen,so(2));
for i=1:so(2)
    [m,n]=min(rsq(:,i));
    stf(:,i)=stfall(:,i,n);
end
return