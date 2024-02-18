function [cx]=gfit0(x,y)
%cx=gfit(x,y),得到x,y的相关系数
s=size(x);
% if s(2)==1
%     a=corrcoef(x,y);
%     aa=min(a);
%     cx=aa(1);
% else
%     for i=1:s(2)
%         a=corrcoef(x(:,i),y);
%         aa=min(a);
%         cx(i)=aa(1);
%     end
% end
nn=sum(x.*y);
mm1=sum(x.*x);
mm2=sum(y.*y);  
if sqrt(mm1.*mm2)==0|nn==0
    cx=0;
    return;
else
    cx=nn./sqrt(mm1.*mm2);
end