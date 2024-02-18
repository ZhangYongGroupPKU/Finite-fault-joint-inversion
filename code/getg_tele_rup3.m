% 断层面参数设置：以下分别设置震中位置（纬度，经度）、断层面解（走向，倾角，滑
% 动角）、子断层划分数（沿倾向，沿走向）、子断层尺寸（沿倾向，沿走向，km），震
% 源所在的子断层位置（沿倾向，沿走向）、green函数所在路径、导入文件名、保存文
% 件名

% get green's function for EVERY SUBFAULT SEPARATELY
function [ob,g,loca,mm,filename] = getg_tele(epi,depth,fault,grid,gridsize,...
    source,gpath,datafolder,leng,srate,range_dis,azim_inter,dist_inter,...
    snr_level,isdis)

% epi = [-42.725,173.065];
% fault = [219,38,128];
% grid = [6,25];
% gridsize = [10,10];
% source = [3,21];
% gpath = 'E:\Studies\joint_invt_pack\green_func\';
% depth = 19;
% datafolder = 'E:\Studies\joint_invt_pack\2016NZ\2016-11-13-mww78-south-island-new-zealand';
% leng = 120;
% srate = 1;
% range_dis = [40,90];
% azim_inter = 5;
% dist_inter = 5;
% snr_level = 10;
% % 若要将速度记录积分为位移记录则为1，否则为0
% isdis = 1;

[obs,mm,loca,chnl,finalpath] = rup_getdata(datafolder,'Z',0,srate,range_dis,epi);
da = distazim_zh(loca,epi);
if max(depth(:))<100
	t = round(get_time(da(:,1)/111.19492664455875,depth,'P')*srate);
else
	t = round(dbgrn_get_time(da(:,1)/111.19492664455875,depth,'P')*srate);
end
nsta = length(t);
ob0 = zeros(650*srate,nsta);
for i = 1:nsta
	obs0 = obs(max(t(i)-50*srate,1):t(i)+600*srate-1,i);
	ob0(1:size(obs0,1),i) = obs0;
end
ob = ob0;
snr = mean(abs(ob0(50*srate+1:100*srate,:))).^2./mean(abs(ob0(1:50*srate,:))).^2;
xx = snr<snr_level;
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];snr(xx) = [];chnl(xx,:) = [];finalpath(xx) = [];
xx = isinf(snr);
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];snr(xx) = [];chnl(xx,:) = [];finalpath(xx) = [];
xx = isnan(snr);
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];snr(xx) = [];chnl(xx,:) = [];finalpath(xx) = [];
ob(1:40*srate,:) = [];

mob = max(abs(ob));
xx = find(mob>mean(mob)*5);
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];chnl(xx,:) = [];finalpath(xx) = [];snr(xx) = [];

[dex] = da_del(da,[dist_inter,azim_inter],'P');
ob = ob(:,dex);loca = loca(dex,:);da = da(dex,:);mm = mm(dex,:);chnl = chnl(dex,:);
finalpath = finalpath(dex);snr = snr(dex);final_snr = snr;

tic
[gt,locasub,dep]=rup3_getg(gpath,fault,grid,gridsize,source,loca,epi,srate);
toc

glen_pre=50; % 100 for sumatra earthquake at teleseismic distance (30 deg, 1200 km long fault), and 150 for near field data inversion
waveleng=700;
nsub=prod(grid);nsta=size(loca,1);
g=zeros(waveleng*srate,6,nsub,nsta);
    
t_begin=max(t-glen_pre*srate,1);
t_end=min(t+waveleng*srate-1-glen_pre*srate,size(gt,1));
for i=1:nsta
    g0=gt(t_begin(i):t_end(i),:,:,i);
    g(1:size(g0,1),:,:,i)=g0;
end

% g has the size of leng*6*nsub*nsta
load([gpath,'\earth.mat']);
miu=earth(:,4).*earth(:,3).^2.*1e9;
[x,y]=histc(dep(:),earth(:,1));
miu=miu(y);
for i=1:nsub
	g(:,:,i,:)=g(:,:,i,:)*miu(i)/3e10;
end

ob(round((leng+glen_pre-10)*srate)+1:end,:) = [];
g(round((leng+glen_pre-10)*srate)+1:end,:,:,:)=[];

if isdis
    g = cumsum(g)/srate;
    ob = cumsum(ob)/srate;
end

filename = 'obg_tele.mat';
save(filename,'ob','g','loca','mm');
end
