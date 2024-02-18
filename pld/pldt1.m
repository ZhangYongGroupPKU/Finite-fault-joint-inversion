function stf=pldt1(a,b,ll,n);
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
bb=b';
b0=[b(1),fliplr(bb(2:end))];
ff=zeros(ll,n);
t=1./max(abs(fft(b,2.^nextpow2(length(b))))).^2;
for i=2:n
    x=conv(ff(:,i-1),b);
    xx=a(1:min(l,length(x)))-x(1:min(l,length(x)));
    xxx=conv(b0',xx);
    ff(:,i)=ff(:,i-1)+t.*xxx(1:ll);
end
stf=ff(:,end);
    