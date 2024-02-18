function [h]=pcolor_zh(x,y,slip,z,flag)
% if nargin==5, the x and y are the coordinates of corners, instead of
% center points

num=numel(slip);

if nargin<5
    x=mat2mat0(x);y=mat2mat0(y);
    x(2:2:end,:)=[];y(2:2:end,:)=[];
    x(:,2:2:end)=[];y(:,2:2:end)=[];
    if nargin>3
        z=mat2mat0(z);
        z(2:2:end,:)=[];z(:,2:2:end)=[];
    end
end

nm=size(x);

xx=zeros(4,num);yy=zeros(size(xx));
if nargin>3;
    zz=zeros(4,num);
end
for i=1:nm(1)-1
    for j=1:nm(2)-1
        xx(:,(j-1)*(nm(1)-1)+i)=[x(i,j);x(i,j+1);x(i+1,j+1);x(i+1,j)];
        yy(:,(j-1)*(nm(1)-1)+i)=[y(i,j);y(i,j+1);y(i+1,j+1);y(i+1,j)];
        if nargin>3;
            zz(:,(j-1)*(nm(1)-1)+i)=[z(i,j);z(i,j+1);z(i+1,j+1);z(i+1,j)];
        end
    end
end

%if nargin==4;zz(dex)=[];end
%numel(slip)

dex=isnan(slip(:));
xx(:,dex)=[];yy(:,dex)=[];slip(dex)=[];
if nargin>3
    zz(:,dex)=[];
end
if nargin==3
%h=patch(xx,yy,slip(:)','facecolor','flat');


h=patch(xx,yy,slip(:)');
elseif nargin>3
    %h=patch(xx,yy,zz,slip(:)','facecolor','flat');
    h=patch(xx,yy,zz,slip(:)');
end
set(h,'LineStyle','none');
set(h,'edgecolor','none');
set(h,'edgecolor',[1,1,1]*0.4,'linewidth',0.2);
box on