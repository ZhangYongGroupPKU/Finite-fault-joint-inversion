function stf=pldt(a,b,ll,n);
%stf=pldf(a,b,ll,n)
%output:
%stf: the source time function
%input:
%a:   seismogram of the main shock
%b:   seismogram of the aftershock
%ll:  the length of the source time function
%n:   the number of iteration
sa=size(a);sb=size(b);
if sa(1)<sa(2) 
    a=a'; 
end
if sb(1)<sb(2) 
    b=b'; 
end
l=length(a);
b=[b;zeros(l-length(b),1)];
m=[b',b'];
mm=fliplr(m);
g=zeros(l,l);
for i=1:l
    g(i,:)=mm(l+1-i:2*l-i);
end

t=1./max(abs(fft(b,2.^nextpow2(length(b))))).^2;
if ll<length(a)
    g=g(:,1:ll);
elseif ll>length(a)
    ll=length(a);
end
gg=g';
ff=zeros(ll,n);

for i=2:n
    xx=ff(:,i-1)+t*gg*(a-g*ff(:,i-1));
    xx(find(xx<0))=0;
    xx(1)=0;
    ff(:,i)=xx;
end
stf=ff(:,end);