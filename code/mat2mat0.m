function m=mat2mat0(mat)
%

%mat=reshape(1:9,[3,3]);
% if nargin==1
%     num=2;
% end
m0=interm(mat,2);
m=zeros(size(mat)*2+1);
m(2:end-1,2:end-1)=m0;

m(1,:)=2*m(2,:)-m(3,:);
m(end,:)=2*m(end-1,:)-m(end-2,:);

m(:,1)=2*m(:,2)-m(:,3);
m(:,end)=2*m(:,end-1)-m(:,end-2);
return
