function moment=fp2mt_zh(m,x,y,z)
%moment=fp2mt_zh(m,x,y,z)
%m: scalar seismic moment
%x: strike
%y: dip
%z: rake
%moment: seismic moment
m=m(:)';
x=x(:)';
y=y(:)';
z=z(:)';
moment=zeros(6,length(x));

x=x.*pi./180;
y=y.*pi./180;
z=z.*pi./180;

sinx=sin(x);cosx=cos(x);sin2x=2.*sinx.*cosx;cos2x=cosx.*cosx-sinx.*sinx;
siny=sin(y);cosy=cos(y);sin2y=2.*siny.*cosy;cos2y=cosy.*cosy-siny.*siny;
sinz=sin(z);cosz=cos(z);
moment(1,:)=-m.*(cosz.*siny.*sin2x+sinz.*sin2y.*(sinx.*sinx));
moment(2,:)=m.*(cosz.*siny.*cos2x+sinz.*sin2y.*sin2x./2);
%moment(2,1)=moment(1,2);
moment(3,:)=-m.*(cosz.*cosy.*cosx+sinz.*cos2y.*sinx);
%moment(3,1)=moment(1,3);
moment(4,:)=m.*(cosz.*siny.*sin2x-sinz.*sin2y.*cosx.*cosx);
moment(5,:)=m.*(-cosz.*cosy.*sinx+sinz.*cos2y.*cosx);
%moment(3,2)=moment(2,3);
moment(6,:)=m.*sinz.*sin2y;
return
% ================================end======================================
% moment(1,:)=-m.*(cos(z).*sin(y).*sin(2.*x)+sin(z).*sin(2.*y).*(sin(x)).^2);
% moment(2,:)=m.*(cos(z).*sin(y).*cos(2.*x)+(1./2).*sin(z).*sin(2.*y).*sin(2.*x));
% %moment(2,1)=moment(1,2);
% moment(3,:)=-m.*(cos(z).*cos(y).*cos(x)+sin(z).*cos(2.*y).*sin(x));
% %moment(3,1)=moment(1,3);
% moment(4,:)=m.*(cos(z).*sin(y).*sin(2.*x)-sin(z).*sin(2.*y).*(cos(x)).^2);
% moment(5,:)=m.*(-cos(z).*cos(y).*sin(x)+sin(z).*cos(2.*y).*cos(x));
% %moment(3,2)=moment(2,3);
% moment(6,:)=m.*sin(z).*sin(2.*y);