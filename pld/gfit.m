function [cx,cxx]=gfit(x,y)
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
if s(2)==1
    if sqrt(mm1.*mm2)==0
        cxx=0.0001;
    else
        cxx=nn./sqrt(mm1.*mm2);
    end
    %cx=cxx.*(1-sqrt(sum((x-y).^2)./sum((abs(x)+abs(y)).^2)));
    cx=(1+cxx)./2.*(1-sqrt(sum((x-y).^2)./sum((abs(x)+abs(y)).^2)));
    %cx=1-(1-(1+cxx)./2).*(sqrt(sum((x-y).^2)./sum((abs(x)+abs(y)).^2)));
else
    if sqrt(mm1.*mm2)==0
        cxx=0.0001*ones(1,s(2));
    else
        cxx=nn./sqrt(mm1.*mm2);
    end
    cx=(1+cxx)./2.*(1-sqrt(sum((x-y).^2)./sum((abs(x)+abs(y)).^2)));
end