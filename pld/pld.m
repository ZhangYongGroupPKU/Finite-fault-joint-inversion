function [stf,l0,ll,c]=pld(a,b,ll,n,t)
%[stf,l0,cc]=pld(a,b,ll,n)
%output:
%stf: the source time function
%l0:  the best length of the stf
%cc:  the misfit
%input:
%a:   seismogram of the main shock
%b:   seismogram of the aftershock
%ll:  [ll(1),ll(2)]:the length range of the source time function
%n:   the number of iteration

sa=size(a);

c=zeros(length(ll),sa(2));
len=sa(1);
j=1;
%s=zeros(2^nextpow2(sa(1)),length(ll));

fftlen=2^nextpow2(sa(1));
pb=fft(b,fftlen);
for i=1:length(ll)%ll(1:end)
    lenlen=ll(i);
    s=pldf(a,b,lenlen,n);

    ps=fft(s);
    ss=real(ifft(pb.*ps));
    %     zx=a(1:len)-ss(1:len);
    %     xx1=zx'*zx;
    [xx1]=gfit(a(1:len,:),ss(1:len,:));

    c(j,:)=1-xx1;
    j=j+1;
end

cc=c;

for i=1:sa(2)
    c=cc(:,i);
    c0=(c-min(c))./(max(c)-min(c));
    lx=ll./ll(end);
    xz=c0.^2+lx'.^2;
    [x,xx]=min(xz);
    lnn=min(xx,length(c));% '+2'
    stf(:,i)=pldf(a(:,i),b(:,i),ll(lnn),n);
    l0(i)=ll(lnn);
    co(i)=c(lnn); %%  '+2'
%     if nargin>4
%         plot(ll(1:end),c,'k')
%         hold on
%         plot(l0,cc(:,i),'ro','markersize',12,'linewidth',2)
%     end
end
% if nargout==4
%     moment=sum(s)';
% end