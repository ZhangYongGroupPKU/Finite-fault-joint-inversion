function plt_area(latlim,longlim,A0,col1,col2)
%
if nargin<5
    col2='w';
    if nargin<4
        col1=[0.95,0.87,0.73];% or ones(1,3)*0.95 for gray color scales, and [0.9,0.9,0.5] for color map
        if nargin<3
             A0=1e1;
        end
    end
end

load coast_map;
%load m_coasts_zh.mat;

% longitude of the coast
long=ncst(:,1);
% latitude of the coast
lat=ncst(:,2);
%select the area within latlim and longlim:
dex=long>longlim(1)&long<longlim(2)&lat>latlim(1)&lat<latlim(2);
loca=ncst(dex,:);
z=find(dex==1);
[~,area_dex]=histc(z,k);
y=unique(area_dex);

%figure
hold on
for i=1:length(y);
    if Area(y(i))<A0%1e1 land
        continue;
    end
    patch(ncst(k(y(i))+1:k(y(i)+1)-1,1),ncst(k(y(i))+1:k(y(i)+1)-1,2),col1,'linewidth',0.25);
end
for i=1:length(y);
    if Area(y(i))>-A0%0 lakes
        continue;
    end
    patch(ncst(k(y(i))+1:k(y(i)+1)-1,1),ncst(k(y(i))+1:k(y(i)+1)-1,2),col2,'linewidth',0.25);
end
% longlim
% latlim
% cosd(mean(latlim))
set(gca,'xlim',longlim,'ylim',latlim,'dataaspectratio',[1,cosd(mean(latlim)),1],...
    'layer','top','box','on');
% set(gca,'xlim',longlim,'ylim',latlim,...
%     'layer','top','box','on');
return


%--------------------------------------------------------------------------
% read the coast data:
% m_proj('Equidistant Cylindrical','long',[-180 180],'lat',[-90 90]);
% [ncst,Area,k]=mu_coast('h',[path,'\gshhs_i.b']);
% path is the position of the file: gshhs_i.b
%--------------------------------------------------------------------------


% [x0,y0]=m_ll2xy(loca1(:,2),loca1(:,1));
% hold on
% %plot(x0,y0,'k','linewidth',2)
% plot(x0,y0,'k.','markersize',15)