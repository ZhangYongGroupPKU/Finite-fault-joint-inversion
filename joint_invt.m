% This is the main script to conduct the rupture inversion process, which
% allows to adopt data of TELESEISMIC record, STRONG-MOTION record, GPS 
% record and INSAR record JOINTLY or SINGLY. 

% Programmed by Yong Zhang, packed by Pingchuan Wang
% 2018/11 at Peking University
% Updated by Sibo Hua,2021.03.15 at Peking University
% Updated by Sibo Hua,2021.09.06 at Peking University
% Updated by Sibo Hua,2024.01.03 at Peking University
%% Parameters of thefault plane and other information (1)

% Location of the epicenter, epi = [latitude,longitude] (1.1)
close all;
clear all;
epi = [-42.725,173.065];
% Depth of the epicenter in km (1.2)
depth = 19;
% Parameters of the fault plane setting, fault = [strike,dip,rake], all in 
% degree (1.3)
fault = [219,38,128];

% Numbers of subfaults along both direction, grid = numbers of subfaults 
% [along dip,along strike] (1.4)
grid = [8,25];

% Size of each subfault in km, gridsize = length of each subfault [along 
% dip,along strike] (1.5)
gridsize = [10,10];

% Which subfault the source locates, source = the ordinal number [along 
% dip,along strike] of the subfault (1.6)
source = [3,21];

% In fact, the upper edge of the fault plane in our model is on the
% ground, so the ordinal number along dip where the source locates is
% determined by the depth and dip. DO NOT CHANGE UNLESS NECESSARY. (1.7)
%source(1) = round(depth/sind(fault(2))/gridsize(1)+1/2);

% Determine to invert in velocity (isdis = 0) or displacement (isdis = 1)
% form (1.8)
isdis = 1;

% Name of the file the inversion results to save in (1.9)
savefile = 'result_joint.mat';

%% Parameters of types of records (2)
% PARAMETERS OF TELESEISMIC RECORDS
% Weight of the records in the joint inversion. If the records are not
% included, then set 0.(2.1)
weight_tele = 1;
teletp=0;%%%远震格林函数的计算方式：每个子断层到每个台站精准计算（teletp=1）或用震源处的格林函数做路径近似（teletp=0）。
ifastf=0;%%%视震源时间函数反演
weightastf=0.1;
% Bandpass of the records in Hz (2.2)
band_tele = [0.01,0.2];

% Name of file the observed record saved in (2.3)
file_tele = 'ob_tele.mat';

% Path of the Green's function for the inversion (2.4)
gpath_tele = 'D:\matlab2021\codes\jtinew2022.03 (3)\example-2016NZ\Green_tele_NZ\green_func\';

% Allow (set 1) or not allow (set 0) the observed records to shift to a 
% time when they fit the synthetic waves the best after an initial
% inversion, and then conduct a second inversion calculation (2.5)
timeshift_tele = 1;

% Calculate (set 1) the Green's function of the records or not (set 0) if
% the calculation has been done and no need to change (2.6) 
ifcalg_tele = 1;

% PARAMATERS OF STRONG-MOTION RECORDS
% Weight of the records in the joint inversion. If the records are not
% included, then set 0.(2.7)
weight_sm = 1;

% Bandpass of the records in Hz (2.8)
band_sm = [0.02,0.20];

% Name of file the observed record saved in (2.9)
file_sm = 'sm_nz2016.mat';

% Path of the Green's function for the inversion (2.10)
gpath_sm = 'D:\matlab2021\codes\jtinew2022.03 (3)\example-2016NZ\Green_sm_NZ\green_func\';

% Allow (set 1) or not allow (set 0) the observed records to shift to a 
% time when they fit the synthetic waves the best after an initial
% inversion, and then conduct a second inversion calculation (2.11)
timeshift_sm = 1;

% Calculate (set 1) the Green's function of the records or not (set 0) if
% the calculation has been done and no need to change (2.12) 
ifcalg_sm = 1;

% PARAMETERS OF GPS RECORDS
% Weight of the records in the joint inversion. If the records are not
% included, then set 0.(2.13)
weight_gps = 0;

% Name of file the observed record saved in (2.14)
file_gps = 'gps.mat';

% PARAMETERS OF INSAR RECORDS
% Weight of the records in the joint inversion. If the records are not
% included, then set 0.(2.15)
weight_insar = 0;

% Name of file the observed record saved in (2.16)
file_insar = 'insar.mat';

%% Parameters of the inversion calculation (3)

% Total length of source time function of the subfaults in s (3.1)
lensub = 110;

% Maximum rupture time of the subfaults in s (3.2)
lenrup = 30;

% Maximum rupture velocity in km/s (3.3)
vrup = 3;

% Spatio-smooth factor (3.4)
lambda1 = 50;

% Temporal-smooth factor (3.5)
lambda2 = 100;

% Minimum scalar moment restriction factor (3.6)
lambda3 = 0.01;

% Allow (set 1) the rake vary in a ±45 degree range or not (set 0) (3.7)
rakevar =1;%%if use astf inversion，rakevar must be 0;

% Set slip on any edge of the fault plane to 0. 0 for setting no edge to 0,
% 1 for setting all edges except for the upper edge to 0, 2 for setting all
% edges to 0. (3.8)
set_edge0 = 1;

%% Parameters of some of the figures (4)

% Range of the map in Figure 2, or just keep them empty for a default
% setup (4.1)
% lon_smsta = [171,176];
% lat_smsta = [-44,-40];
lon_smsta = [];
lat_smsta = [];

% rct = numbers of [rows,column] in Figure 5, or just keep it empty for a
% default setup (4.2)
rct = [];

% rct = numbers of [rows,column] in Figure 6, or just keep it empty for a
% default setup (4.3)
rcsm = [];

% Range of the map in Figure 7, or just keep them empty for a default
% setup (4.4)
lon_gps = [172,177];
lat_gps = [-44,-40];

% Ticks interval (in degree) of the maps in Figure 8 (4.5)
dlon_sar = 0.4;
dlat_sar = 0.4;

%% Inversion calculations (5)
[miu1,miun] = getmiu(depth,fault,grid,gridsize,source,epi);%断层上各处的剪切模量？
if ifastf==1
    [obt,gt,locat,mmt,s1,leng,t,t0] = getobg_tele4(weight_tele,ifcalg_tele,file_tele,...
    epi,depth,fault,grid,gridsize,source,gpath_tele,miun,band_tele,isdis);
    tdt=t-repmat(t0',size(t,1),1);
    
    
    if teletp==1
[obt,gt,locat,mmt,s1,leng] = getobg_tele2(weight_tele,ifcalg_tele,file_tele,...
     epi,depth,fault,grid,gridsize,source,gpath_tele,miun,band_tele,isdis);%%%%%精准子断层&台站格林函数
    end
elseif teletp==0
 [obt,gt,locat,mmt,s1,leng] = getobg_tele3(weight_tele,ifcalg_tele,file_tele,...
     epi,depth,fault,grid,gridsize,source,gpath_tele,miun,band_tele,isdis);%%%%考虑路径近似的格林函数
elseif teletp==1
    [obt,gt,locat,mmt,s1,leng] = getobg_tele2(weight_tele,ifcalg_tele,file_tele,...
     epi,depth,fault,grid,gridsize,source,gpath_tele,miun,band_tele,isdis);
    
    
end


[obsm,gsm,locasm,mmsm,s2] = getobg_sm(weight_sm,ifcalg_sm,file_sm,...
    epi,depth,fault,grid,gridsize,source,gpath_sm,miun,band_sm,isdis);%强震观测记录、格林函数、台站经纬度、名称、采样率

[obt,gt] = length_unification(obt,gt,obsm);%有效记录长度统一，短的记录补0

srate = set_srate(s1,s2);%验证采样率统一

ob = [obt,obsm];g = cat(4,gt,gsm);loca = [locat;locasm];mm = [mmt;mmsm];%连接矩阵
clear gt gsm

[ob,g,mob] = normalize_obg(ob,g,obt,obsm,weight_tele,weight_sm);%平权和加权

[G,sGDT0,hdata,meangreen,sumob,trup_con,locasub,locadep,trup] = get_Gmatrix(...
    ob,g,loca,epi,depth,fault,grid,gridsize,source,lensub,lenrup,vrup,srate,...
    lambda1,lambda2,lambda3,rakevar,set_edge0);%反演有限断层矩阵构建
if ifastf==1
    obtele=ob(:,1:size(locat,1));
    gtele=g(:,:,:,1:size(locat,1));
    [astftele,synastf]=getastftele(obtele,gtele,source,grid,fault,lensub,band_tele,s1);
    
    
    Gastf=get_stfinvG(grid,trup_con,tdt,lensub*srate);
    astftele=astftele(:)/srate;
    mobstf=max(max(Gastf))/weightastf/max(max(G))*2;
    astftele=astftele/mobstf;
    Gastf=Gastf/mobstf;
    grida=grid;
    save('astfneed.mat','gtele','source','grida','fault','lensub','band_tele','s1','rakevar','mobstf');
else
    Gastf=[];
    astftele=[];
end
[slipg,Ggps,weit_gps,locagps] = getobg_gps2(weight_gps,file_gps,sumob,...
    locasub,locadep,fault,grid,gridsize,trup_con,rakevar);%GPS滑动、格林函数、权重、位置

[slipsar,Gsar,weit_sar,locasar,nsar,dr] = getobg_insar2(weight_insar,file_insar,...
    sumob,locasub,locadep,fault,grid,gridsize,trup_con,rakevar);%insar滑动、格林函数、权重、位置

G_jt = [G;Ggps*weit_gps;Gsar*weit_sar];%反演格林函数矩阵构建
clear g G Ggps Gsar
sGDT=[Gastf;sGDT0];
hdata_use = [ob(:)/meangreen;slipg(:)*weit_gps;slipsar*weit_sar;astftele(:);zeros(...
    size(sGDT0,1),1)];%反演用观测记录矩阵构建
tic;[Xsolve,rsq] = cgls_xu_invt0(G_jt,sGDT,hdata_use,1e-170,1000);toc%反演
syn = reshape(G_jt(1:numel(ob),:)*(Xsolve*meangreen),size(ob));%合成记录
%plot(rsq);
%saveas(gcf,'rsq1','fig');
[Xsolve,ob0,syn,astftele2,synastf2] = shift_ob(timeshift_tele*weight_tele,timeshift_sm*...
    weight_sm,Xsolve,ob,syn,obt,obsm,hdata_use,meangreen,G_jt,sGDT,srate,ifastf);%时移后再次计算
if ifastf==1
synstf=sGDT(1:size(obt,2)*lensub*srate,:)*Xsolve*mobstf*srate;
%synastf2=synastf2;
end

misfta0=sum(sum((ob0-syn).^2))/sum(sum(ob0.^2));

[X,ob,syn] = get_obgX(Xsolve,ob0,syn,mob,lensub,trup_con,prod(grid),srate,rakevar);%%%%try

obsm = ob(:,size(obt,2)+1:end);
synsm = syn(:,size(obt,2)+1:end);
obsm = reshape(obsm,[size(ob,1),size(obsm,2)/3,3]);
obsm = permute(obsm,[1,3,2]);
synsm = reshape(synsm,[size(ob,1),size(obsm,3),3]);
synsm = permute(synsm,[1,3,2]);
syng = G_jt(numel(ob)+1:numel(ob)+numel(slipg),:)/weit_gps*Xsolve;
synsar = G_jt(numel(ob)+numel(slipg)+1:end,:)/weit_sar*Xsolve;

%[X,ob,syn] = get_obgX(Xsolve,ob,syn,mob,lensub,trup_con,prod(grid),srate);

%[slip,stf,slipx,slipy,slipe,slipn,zslip] = get_result(X,lensub,grid,...
%    gridsize,fault,miun,srate);
[slip,stf,slipx,slipy,slipe,slipn,zslip,stfsub_out] = get_result2(X,lensub,grid,...
    gridsize,fault,miun,srate);
stfsub_out=stfsub_out*srate;
mw = m2m(sum(slip(:).*miun(:)*prod(gridsize)*1e6));

msftj=sum(sum((ob-syn).^2))/sum(sum((ob).^2));
if weight_tele
obtelef=ob(:,1:size(obt,2));
syntele=syn(:,1:size(obt,2));
msftt=sum(sum((obtelef-syntele).^2))/sum(sum((obtelef).^2));
msftti=sum((obtelef-syntele).^2)./sum((obtelef).^2);
end
if weight_sm
msftsm=sum(sum(sum((obsm-synsm).^2)))/sum(sum(sum((obsm).^2)));
msftsmi(:,:)=sum((obsm-synsm).^2)./sum((obsm).^2);
end
%% Figures (6)

clear G_jt G g ncst sGDT k Area Gsar G_gps
save(savefile)
[dex,dey]=substfs_plot(epi,grid,gridsize,source,stfsub_out,srate);
%saveas(gcf,'stfsub','fig');
% Figure 2: Distribution of the teleseismic stations
if weight_tele
    figure
    epi_sta00(epi,loca(1:size(obt,2),:),mm(1:size(obt,2),:),100);
    title('Teleseismic Stations');
    %saveas(gcf,'locatele','fig');
end

% Figure 3: Distribution of the strong-motion stations
if weight_sm
    plt_sm_sta(lon_smsta,lat_smsta,epi,locasm,mmsm);
    %saveas(gcf,'locasm','fig');
end

% Figure 4: Static slip distribution
figure
[numth,interval,m] = colornumth(slip);
[slip01,start,stop] = faultslip_pcolor(grid,gridsize,source,zslip,...
    (numth)*interval,numth,interval);
axis equal tight
set(gca,'linewidth',0.5,'tickdir','in','xminortick','off','yminortick',...
    'off','fontsize',16);
%saveas(gcf,'rup','fig');

% Figure 5: Source time function
if ~isempty(ob)
    figure
    fill(0:1/srate:lensub,[stf;0],'c');
    set(gca,'layer','top');
    set(gca,'fontsize',15);
    title(['Source Time Function',num2str(mw,'%.2f')]);
    xlabel('Time(s)');
    ylabel('Moment rate(Nm/s)');
    %saveas(gcf,'stf','fig');
end

% Figure 6: Fitting of teleseismic waves
if weight_tele
    if isempty(rct)
        rct = [4,ceil(size(obt,2)/4)];
    end
    leng=leng*srate;
    pltobsyn_new(ob(1:leng,1:size(obt,2)),syn(1:leng,1:size(obt,2)),...
        (1:leng)/srate,rct,0.05,mm(1:size(obt,2),:),repmat('P/Z',...
        [size(obt,2),1]));
    %saveas(gcf,'syntele','fig');
end

% Figure 7: Fitting of strong-motion waves
if weight_sm
    if isempty(rcsm)
        rcsm = [4,ceil(size(obsm,3)/4)];
    end
    pltobsyn3_new(obsm,synsm,0:1/srate:(size(obsm,1)-1)/srate,rcsm,0.05,...
        mm(size(obt,2)+1:size(obt,2)+size(obsm,3),:));
    %saveas(gcf,'synsm','fig');
end

% Figure 8: Fitting of GPS records
if weight_gps
    matching_gps(lon_gps,lat_gps,locagps,locasub,epi,grid,slip,slipg,syng);
    title('GPS Data Fitting');
end

% Figure 9: Fitting of InSAR records
if weight_insar
%     plot([slipsar,synsar])
    matching_insar(dlon_sar,dlat_sar,slipsar,synsar,locasar,nsar,dr);
%     title('InSAR Data Fitting');
end

disp(['The Mw calculated is ',num2str(mw,'%.2f'),'!'])

n1=8;%时间间隔（s）
n2=3;%图行数
n3=5;%图列数
stfsub_temp=stfsub_out./miun(:)'*3e10;
stfsub_sliprate=stfsub_temp(:)/3e16/prod(gridsize);
maxslip=max(max(slip));
%Figure 10 每n1秒时间段内的累计滑动量
rate=rup_pltslipshot2(stfsub_sliprate/srate,lensub*srate,[n1*srate,n2,n3],1,maxslip/2,12,srate);
colorbar('southoutside','position',[0.27    0.0516    0.5    0.0308]);
%Figure 11 截止到不同个n1秒时间的累计滑动量
rate=rup_pltslipshot2(stfsub_sliprate/srate,grid,lensub*srate,[n1*srate,n2,n3],2,maxslip,12,srate);
colorbar('southoutside','position',[0.27    0.0516    0.5    0.0308]);
%%Figure 12、13 视震源时间函数反演拟合
if ifastf==1
    if isempty(rct)
        rct = [4,ceil(size(obt,2)/4)];
    end
    %leng=leng*srate;
    if size(synastf2,1)>0
        synastf=synastf2;
        astftele=astftele2;
    end
    astftele=reshape(astftele,lensub*srate,size(obt,2));
    synstf=reshape(synstf,lensub*srate,size(obt,2));
    astftele=astftele*mobstf*srate;
    
    pltobsyn_new(ob0(1:leng,1:size(obt,2)),synastf(1:leng,1:size(obt,2)),...
        (1:leng)/srate,rct,0.05,mm(1:size(obt,2),:),repmat('P/Z',...
        [size(obt,2),1]));
    rct = [4,ceil(size(obt,2)/4)];
    pltobsyn_new(astftele(1:lensub*srate,1:size(obt,2)),synstf(1:lensub*srate,1:size(obt,2)),...
    (1:lensub*srate)/srate,rct,0.05,mm(1:size(obt,2),:),repmat('P/Z',...
        [size(obt,2),1]))
end
