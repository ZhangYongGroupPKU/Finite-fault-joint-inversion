% This is the script to pick, repick and view teleseismic waves

%% 
% Depth of the epicenter in km
depth = 19;

% Length of the wave before initial P in s
leng_before = 30;

% Length of the wave after initial P in s
leng_after = 120;

% Sample rate in s^-1
srate = 1;

% Location of the epicenter, epi = [latitude,longitude]
epi = [-42.725,173.065];

% Range of epicentral distance (in degree) within which the stations could
% be selected, as teleseismic distances 
range_distance = [40,90];

% Minimum take-off angle interval (in degree) between each stations could be 
% selected
azim_inter = 5;

% Minimum epicentral distance interval (in degree) between each stations
% could be selected
dist_inter = 5;

% Minimum signal-to-noise-ratio (in %) under which the wave would not be 
% selected
snr_level = 10;

% Path of teleseismic records downloaded from IRIS 
datafolder = 'D:\matlab2021\codes\jtinew2022.03 (3)\example-2016NZ\2016-11-13-mww78-south-island-new-zealand';

% Ordinal numbers of stations that are to exclude from the inversion
deldex = [];

% Decide whether to show filtered (set 1) or original (set 0) waves in the
% following figure
iffilter = 1;

% Bandpass (in Hz) of the records shown in the following figure
% Bandpass here WILL NOT BE SAVED AND USED in the inversion
band = [0.01,0.2];

%%
leng = leng_before + leng_after;

[obs,mm,loca,chnl,finalpath] = rup_getdata(datafolder,'Z',0,srate,range_distance,epi);

da = distazim_zh(loca,epi);
if max(depth(:))<100
	t = round(get_time(da(:,1)/111.19492664455875,depth,'P')*srate);
else
	t = round(dbgrn_get_time(da(:,1)/111.19492664455875,depth,'P')*srate);
end
nsta = length(t);
ob0 = zeros(650*srate,nsta);
for i = 1:nsta
	obs0 = obs(max(t(i)-leng_before*srate,1):t(i)+(650-leng_before)*srate-1,i);
	ob0(1:size(obs0,1),i) = obs0;
end
ob = ob0;
snr = mean(abs(ob0(leng_before*srate+1:(leng_before+50)*srate,:))).^2./mean(abs(ob0(1:leng_before*srate,:))).^2;
xx = snr<snr_level;
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];snr(xx) = [];chnl(xx,:) = [];finalpath(xx) = [];
xx = isinf(snr);
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];snr(xx) = [];chnl(xx,:) = [];finalpath(xx) = [];
xx = isnan(snr);
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];snr(xx) = [];chnl(xx,:) = [];finalpath(xx) = [];
% ob(1:40*srate,:) = [];

mob = max(abs(ob));
xx = find(mob>mean(mob)*5);
ob(:,xx) = [];loca(xx,:) = [];da(xx,:) = [];mm(xx,:) = [];chnl(xx,:) = [];finalpath(xx) = [];snr(xx) = [];

[dex] = da_del(da,[dist_inter,azim_inter],'P');
ob = ob(:,dex);loca = loca(dex,:);da = da(dex,:);mm = mm(dex,:);chnl = chnl(dex,:);
finalpath = finalpath(dex);snr = snr(dex);final_snr = snr;

ob(leng*srate+1:end,:) = [];

ob(:,deldex) = [];loca(deldex,:) = [];da(deldex,:) = [];mm(deldex,:) = [];
snr(deldex) = [];chnl(deldex,:) = [];finalpath(deldex) = [];

obf = ob;
if iffilter
    [bb,aa] = butter(3,band*2/srate);
    obf = filter(bb,aa,ob);
end

nsta = size(ob,2);
disp([num2str(nsta),' stations have been reselected!']);

figure
h = epi_sta00(epi,loca,mm,100);

figure
xx = 0:1/srate:length(obf(:,1))/srate;
ss = size(obf);
astf = obf;
astf = [astf;zeros(1,ss(2))];
for i = 1:ss(2)
    astf(:,i) = astf(:,i)./max(abs(astf(:,i)))+i*2;
end
nplt = ceil(nsta/10);
mda = num2str(round([da(:,1)/6371/pi*180,da(:,2)]));
mda = [repmat('(',ss(2),1),mda(:,1:2),repmat('/',ss(2),1),mda(:,4:end),repmat(')',ss(2),1)];
zzchar = [mm,mda];
for i = 1:size(zzchar,1)
    zz0 = zzchar(i,:);
    zz0(zz0 == ' ') = [];
    yl{i} = zz0;
end
for i = 1:nplt
    subplot(1,nplt,i)
    plot(xx,astf(:,(i-1)*10+1:min([i*10,nsta])),'color',ones(1,3)*0.0,'linewidth',1)
    axis tight
    xlabel('Time(s)')
    grid on
    set(gca,'ytick',2:2:ss(2)*2);
    set(gca,'yticklabel',yl,'fontsize',8);
end

save ob_tele.mat ob loca da mm srate leng leng_before
