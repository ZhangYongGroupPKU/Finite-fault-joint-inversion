function [dt,m,rm]=bestco_zh(x,y)
%
k=max(length(x),length(y));
Rxy=xcorr(x,y);
Rxx0=x'*x;
Ryy0=y'*y;
rm=Rxy/sqrt(Rxx0*Ryy0);
%m=(-k+1):(k-1);
[mx,nx]=max((rm));
%dt=m(nx);
dt=nx-k;
m=mx;
