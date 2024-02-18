function plt_sm_sta(lon_smsta,lat_smsta,epi,locasm,mmsm)
if isempty(lon_smsta) || isempty(lat_smsta)
    plt_lon = [min(locasm(:,2)),max(locasm(:,2))] + [-0.5,0.5];
    plt_lat = [min(locasm(:,1)),max(locasm(:,1))] + [-0.5,0.5];
else
    plt_lon = lon_smsta;
    plt_lat = lat_smsta;
end

figure
plt_area(plt_lat,plt_lon);

lontick = ceil(plt_lon(1)):1:floor(plt_lon(2));
lattick = ceil(plt_lat(1)):1:floor(plt_lat(2));

for i = 1:size(mmsm,1)/3
    text(locasm(i,2),locasm(i,1),['   ',mmsm(i,:)],'fontsize',12);
    plot(locasm(i,2),locasm(i,1),'k^','markersize',8,'markerface','c');
    hold on
end
plot(epi(2),epi(1),'kp','markersize',20,'markerface','y')
set(gca,'dataaspectratio',[1,cosd(epi(1)),1]);

[lonlabel,latlabel] = tick2label(lontick,lattick);
set(gca,'xtick',lontick);set(gca,'xticklabel',lonlabel);
set(gca,'ytick',lattick);set(gca,'yticklabel',latlabel);
set(gca,'fontsize',15);
end
