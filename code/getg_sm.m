% 断层面参数设置：以下分别设置震中位置（纬度，经度）、断层面解（走向，倾角，滑
% 动角）、子断层划分数（沿倾向，沿走向）、子断层尺寸（沿倾向，沿走向，km），震
% 源所在的子断层位置（沿倾向，沿走向）、green函数所在路径、导入文件名、保存文
% 件名
function [ob,g,loca,mm,filename] = getg_sm(epi,fault,grid,gridsize,...
    source,gpath,file_sm,srate,isdis)

% epi = [-42.725,173.065];
% fault = [219,38,128];
% grid = [6,25];
% gridsize = [10,10];
% source = [3,21];
% gpath = 'E:\Studies\2016NZ_joint_for_Wang\gdb\green_func\';
% file_sm = 'sm_nz2016.mat';
% isdis = 1;
% srate = 1;

load(file_sm);
leng = size(ob,1);
[locasub,dep] = get_subloca(fault,grid([2,1]),gridsize([2,1]),[1,source([2,1])],loca,epi);
%return

nsta = size(loca,1);
nsub = size(locasub,1);
locasuba = repmat(locasub,[nsta,1]);
depa = repmat(dep,[nsta,1]);
locaa = repmat_zh(loca,nsub);% [... nsub,nsta]

daa = da_zh(locaa,locasuba,2);

dda = [depa,daa];
tic
M = fp2mt_zh(1,fault(1),fault(2),fault(3));
bda = da_zh(locasuba,locaa,2);% get the back azimuth of the subfaults relative to stations
bfai = bda(:,2)+180; % the ray angle
sfai = sin(bfai*pi/180);cfai = cos(bfai*pi/180);
for i = 6:-1:1
    M = zeros(6,1);M(i) = 1;
    g0 = seekg_wang(gpath,dda,M,srate);
    for j = 1:size(g0,3)
        v = g0(:,1,j);r = g0(:,2,j);t = g0(:,3,j); %
        g0(:,1,j) = r*sfai(j)-t*cfai(j); % E
        g0(:,2,j) = r*cfai(j)+t*sfai(j); % N
        g0(:,3,j) = v; % U
    end
    
    g1 = reshape(g0(:,1,:),[size(g0,1),nsub,nsta]);
    g2 = reshape(g0(:,2,:),[size(g0,1),nsub,nsta]);
    g3 = reshape(g0(:,3,:),[size(g0,1),nsub,nsta]);
    
    g(:,i,:,:) = cat(3,g1,g2,g3); % g has the size of [leng,6,nsub,nsta]
end
toc
clear g0 g1 g2 g3

% g has the size of leng*6*nsub*nsta
load([gpath,'\earth.mat']);
miu = earth(:,4).*earth(:,3).^2.*1e9;
[x,y] = histc(dep(:),earth(:,1));
miu = miu(y);
for i = 1:nsub
    g(:,:,i,:) = g(:,:,i,:)*miu(i)/3e10;
end

ob = permute(ob,[1,3,2]);
ob = ob(:,:);
ob(leng+1:end,:) = [];
g(leng+1:end,:,:,:) = [];

if isdis == 1;
    g = cumsum(g)/srate;
%     ob = cumsum(ob)/srate;
end

for i = 1:isdis+1
    ob = cumsum(ob)/srate;
end

filename = 'obg_sm.mat';
save(filename,'ob','g','loca','mm');
