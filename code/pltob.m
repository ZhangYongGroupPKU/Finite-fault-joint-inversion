function pltob(ob,rc,mm,srate)
%==========================================================================
% ob: [leng,3,nsta]
% syn: [leng,3,nsta]
% t: 0:1/srate:size(ob,1)/srate-1
% rc: 
%==========================================================================
% rc
% rat
if length(size(ob)) == 2
    ob = reshape(ob,length(ob),3,[]);
end
    
t = 0:1/srate:(size(ob,1)-1)/srate;
rat = 0.05;

wid0=1/(rc(2)*(1+rat)+rat);
wid=wid0*rat;

fsize=7;
rat_sta=1.0;
nsta=size(ob,3);
rc1=rc(1);
ylim=[rc1*6+(rc1-1)*rat_sta,0];
yzero=zeros(3,rc1);
for j=1:rc1
    yzero(:,j)=ylim(1)-[1;3;5]-(j-1)*rat_sta-6*j;
end
   yzero=yzero(:);
    
figure('color','w')

nt=min(t);mt=max(t);
for i=1:rc(2)
    axes('position',[(i-1)*wid0+wid*i,0.08,wid0,0.9]);
    
    dex0=(i-1)*rc1+1;
    if dex0>nsta;
        axis off
        return;
    end
    
    ob0=ob(:,:,(i-1)*rc1+1:min(i*rc1,end));
    mm0=mm((i-1)*rc1+1:min(i*rc1,end),:);
    
    [ob0,~]=vecnor_w(ob0(:,:),1);
    for j=1:size(ob0,2)
        ob0(:,j)=ob0(:,j)+yzero(j);
    end    
    
    plot(t,ob0,'k','linewidth',0.8);%1.0
    hold on
     
    for j=1:size(mm0,1)
        text(nt,yzero((j-1)*3+1),mm0(j,:),'fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','bottom');
        if i==1&&j==1
            text(nt,yzero((j-1)*3+1),'EW','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
            text(nt,yzero((j-1)*3+1)-2,'NS','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
            text(nt,yzero((j-1)*3+1)-4,'UD','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
        end
    end
    axis tight
    set(gca,'ycolor','w','ytick',[],'tickdir','out','xminortick','on',...
        'ylim',[-7,max(yzero)]+1,'ticklength',[0.02,0.02]);%,'grid','on')
    box off
    xlabel('Time(s)')
end

%==========================================================================
function [mout1,mm]=vecnor_w(min1,rat)
%
if nargin==2;
    rat=1;
end
sm=size(min1);
mm1=max(abs(min1));
mm=mm1;
for i=1:sm(2)
    if mm(i)==0
        continue;
    end
    min1(:,i)=min1(:,i)/mm(i)*rat;%-i*2;
end
mout1=min1;