function [cxy]=gfit1(x,y)
%cxy=gfit1(x,y),得到x,y的相关系数
%cxy is the correlation coeficient of x and y

s=size(x);
if s(2)==1
    xy=x'*y;
    xx=x'*x;
    yy=y'*y;
    if xx==0||yy==0
        cxy=0;
        return;
    else
        cxy=xy./sqrt(xx.*yy);
    end
    return
else
    xy=sum(x.*y);
    xx=sum(x.*x);
    yy=sum(y.*y);
    
    cxy=zeros(size(xy));
    sxy=sqrt(xx.*yy);
    cxy(sxy==0)=0;
    
    nz=sxy~=0;
    cxy(nz)=xy(nz)./sxy(nz);
end