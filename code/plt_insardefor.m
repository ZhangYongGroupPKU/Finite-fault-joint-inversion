% locasar=insar(:,1:2);
% slipsar=insar(:,3);
function plt_insardefor(loca,deform,mndeform,dr,clb)

%clear all
%close all

% col=jet(1e3);
% dex=round((deform-mndeform(2))/(mndeform(1)-mndeform(2))*999)+1;
% hold on
% for i=1:size(loca,1)
%     plot(loca(i,2),loca(i,1),'ko','markerfacecolor',col(dex(i),:));
% end
if nargin<4
numdata=size(loca,1);
loca1=repmat(loca,numdata,1);
loca2=repmat_zh(loca,numdata);
da=da_zh(loca1,loca2);
dist=da(:,1);
dist(1:numdata+1:end)=inf;
dr=min(dist)/2/111.1949;
end

%xx=zeros(4,num);yy=zeros(size(xx));
[A1,A2]=reckon(loca(:,1),loca(:,2),dr*sqrt(2),315);
[B1,B2]=reckon(loca(:,1),loca(:,2),dr*sqrt(2),45);
[C1,C2]=reckon(loca(:,1),loca(:,2),dr*sqrt(2),135);
[D1,D2]=reckon(loca(:,1),loca(:,2),dr*sqrt(2),225);
xx=[A2,B2,C2,D2]';
yy=[A1,B1,C1,D1]';

% size(xx)
% size(yy)
% size(deform)
zp=patch(xx,yy,deform');
set(zp,'edgecolor','w');set(zp,'edgecolor','none')
caxis([mndeform]);
set(gca,'dataaspectratio',[1,cosd(mean(loca(:,1))),1])
colormap(jet_zh2017(100));
axis tight
axis fill
if clb
    h = colorbar;
    set(h,'fontsize',12);
    set(get(h,'Title'),'string','m');
else
end
box on

