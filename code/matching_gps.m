function matching_gps(lon_gps,lat_gps,locagps,locasub,epi,grid,slip,slipg,syng)
if isempty(lon_gps) || isempty(lat_gps)
    plt_lon = [min(locagps(:,2)),max(locagps(:,2))] + [-0.5,0.5];
    plt_lat = [min(locagps(:,1)),max(locagps(:,1))] + [-0.5,0.5];
else
    plt_lon = lon_gps;
    plt_lat = lat_gps;
end

figure
plt_area(plt_lat,plt_lon);

hold on
xx = reshape(locasub(:,2),grid);
yy = reshape(locasub(:,1),grid);
zz = slip;
pcolor_zh(xx,yy,zz);
colormap(jet_zh(100,2))
h = colorbar;
set(get(h,'title'),'string','m');
hold on
quiver(locagps(:,2),locagps(:,1),slipg(:,1),slipg(:,2),1,'k');
syng = reshape(syng,size(slipg));
quiver(locagps(:,2),locagps(:,1),syng(:,1),syng(:,2),1,'r');
load coasts_large.mat
plot(ncst(:,1),ncst(:,2),'k');
plot(epi(2),epi(1),'kp','markersize',20,'markerfacecolor','w');
set(gca,'dataaspectratio',[1,cosd(epi(1)),1]);

lontick = ceil(plt_lon(1)):1:floor(plt_lon(2));
lattick = ceil(plt_lat(1)):1:floor(plt_lat(2));
[lonlabel,latlabel] = tick2label(lontick,lattick);
set(gca,'xtick',lontick);set(gca,'xticklabel',lonlabel);
set(gca,'ytick',lattick);set(gca,'yticklabel',latlabel);
set(gca,'fontsize',15);
end
