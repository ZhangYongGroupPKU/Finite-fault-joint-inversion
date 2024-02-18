function h=epi_sta_zh(epi,sta,m,radius)
%epi_sta([23,101],str2num(mm(:,6:25)),mm(:,1:4),100);
%h=epi_sta(radius)
if nargin==3
    radius=180;
end
%m_proj('stereographic','lat',epi(1),'long',epi(2),'radius',radius);
%m_proj('Azimuthal Equal-area','lat',epi(1),'long',epi(2),'radius',radius);
%figure
m_proj('Azimuthal Equidistant','lat',epi(1),'long',epi(2),'radius',radius);
%m_proj('hammer','clong',epi(2));                                           % projecting methods
m_coast('patch',[.7 .7 .7]*8/7,'edgecolor',[.7 .7 .7]*8/7);                         % color of the continent 
%m_coast('patch',[1 .95 .8],'edgecolor','none'); 
%m_coast('patch',[1 0.95 0.9]*7/7,'edgecolor',[.7 .7 .7]*0/7);                         % color of the continent 
%m_grid('linest','-','xticklabels',[],'yticklabels',[]);
m_grid('xtick',[],'ytick',[],'linestyle','-');                              % circle    
axis equal
axis off
hold on
%m_plot(epi(2),epi(1),'rh','markersize',24,'markerfacecolor','r');           % plot the epicenter
%m_plot(epi(2),epi(1),'r.','markersize',34,'markerfacecolor',[0.3 0.2 0]);           % plot the epicenter
[epi2,epi1]=m_ll2xy(epi(2),epi(1));
%aniseStar([epi2,epi1],0.15,1.1,'w');
pltstar(8,[epi2,epi1],0.15,0.48,'w','k',1);
ss=size(sta);
for i=1:ss(1)
    %m_plot(sta(i,2),sta(i,1),'^','markersize',6,'markerfacecolor',[0.2,0.2,0.2]*0,'markeredgecolor',[0.2,0.2,0.2]*0);   % stations
    m_plot(sta(i,2),sta(i,1),'k^','markersize',10,'markerfacecolor','c');   % stations
    m_text(sta(i,2),sta(i,1),m(i,:),'fontsize',10,'fontname','times new roman');                        % station labels 
end                                                                         %plot stations
%m_range_ring(epi(2),epi(1),[45:45:180]*111.2,'color','k','linewi',0.5);     %plot the distance circle 
m_range_ring(epi(2),epi(1),radius*111.1945,'color','k','linewidth',0.5);     %plot the distance circle 
h=gca;
%set(h,'backgroundcolor','r')
if nargout==0
    clear h;
end
camzoom(1.2)
return
% projecting={'Stereographic',...
%      'Orthographic'... 
%      'Azimuthal Equal-area'... 
%      'Azimuthal Equidistant'... 
%      'Gnomonic'... 
%      'Satellite'... 
%      'Albers Equal-Area Conic'... 
%      'Lambert Conformal Conic'... 
%      'Mercator'... 
%      'Miller Cylindrical'... 
%      'Equidistant Cylindrical'... 
%      'Oblique Mercator'... 
%      'Transverse Mercator'... 
%      'Sinusoidal'... 
%      'Gall-Peters'... 
%      'Hammer-Aitoff'... 
%      'Mollweide'... 
%      'UTM'};