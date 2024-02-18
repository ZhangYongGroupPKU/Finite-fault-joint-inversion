function [crust1] = crust1model(Lat,Lon,sediment)
%基于Crust1.0模型，生成能用于dbgrn的单点速度模型文件
if nargin < 3
    sediment = 0;
end

Lon = Lon + 0.1*(mod(Lon,1)==0);
Lat = Lat + 0.1*(mod(Lat,1)==0);
[~,~,data0] = imp(Lon-0.5,Lon+0.5,Lat-0.5,Lat+0.5);
% data0 = reshape(data0,[9,4]);

data0(1:2,:) = [];
data0(:,1) = data0(1,1) - data0(:,1);

if ~sediment
    data0(1:3,:) = [];
end

data0(data0(:,2)==0,:) = [];

ndep = size(data0,1);
depth = zeros(2*ndep-1,1);
vp = zeros(2*ndep-1,1);
vs = zeros(2*ndep-1,1);
rho = zeros(2*ndep-1,1);
depth(1) = 0;
rho(1) = data0(1,2);
vp(1) = data0(1,3);
vs(1) = data0(1,4);

for i = 2:2*ndep-1
    depth(i) = data0(floor(i/2)+1,1);
    rho(i) = data0(ceil(i/2),2);
    vp(i) = data0(ceil(i/2),3);
    vs(i) = data0(ceil(i/2),4);
end
data = [depth,vp,vs,rho];

crust = load('crust.txt');
crust(:,1) = [];
qp = 927.34*ones(size(data,1),1);
qs = 599.99*ones(size(data,1),1);

crust1 = vertcat([data,qp,qs],crust(6:end,:));

FileName = sprintf('Crust1.0_Lon%.1f_Lat%.1f.txt',Lon,Lat);
fid = fopen(FileName,'wt');
for i = 1:size(crust1,1)
    fprintf(fid,' %g',i);
    fprintf(fid,' %.2f %.4f %.4f %.4f',crust1(i,1),crust1(i,2),...
crust1(i,3),crust1(i,4));
    fprintf(fid,' %.2f %.2f\n',crust1(i,5),crust1(i,6));
end

% create earth.mat
load earth0.mat;
earth = vertcat(data,earth(6:end,:));
save earth.mat earth
end
