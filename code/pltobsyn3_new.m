function [fit] = pltobsyn3_new(ob,syn,t,rc,rat,mm)
%==========================================================================
% ob: [leng,3,nsta]
% syn: [leng,3,nsta]
% t: 0:1/srate:size(ob,1)/srate-1
% rc: 
%==========================================================================
% rc
% rat

wid0=1/(rc(2)*(1+rat)+rat);
wid=wid0*rat;
%high=1/rc(1);

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
    
fit=gfit1(ob,syn);
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
    syn0=syn(:,:,(i-1)*rc1+1:min(i*rc1,end));
    mm0=mm((i-1)*rc1+1:min(i*rc1,end),:);
    
    [ob0,syn0,maxwave]=vecnor(ob0(:,:),syn0(:,:),1);
    for j=1:size(ob0,2)
        ob0(:,j)=ob0(:,j)+yzero(j);
        syn0(:,j)=syn0(:,j)+yzero(j);
    end    
%     ob0=reshape(ob0,[size(ob0,1),3,size(ob0,2)/3]);
%     syn0=reshape(syn0,[size(ob0,1),3,size(ob0,2)/3]);
    
     plot(t,ob0,'k','linewidth',0.8);%1.0
     hold on
     plot(t,syn0,'r','linewidth',0.5);%0.5
    
%     for k=1:size(ob0,2)
%         m3=mod(k,3);
%         if m3==1
%             plot(t,syn0(:,k),'b')
%         elseif m3==2
%             plot(t,syn0(:,k),'r')
%         else
%             plot(t,syn0(:,k),'c')
%         end
%     end
    
    
    for j=1:size(mm0,1)
        text(nt,yzero((j-1)*3+1),mm0(j,:),'fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','bottom');
        if i==1&&j==1
            text(nt,yzero((j-1)*3+1),'EW','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
            text(nt,yzero((j-1)*3+1)-2,'NS','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
            text(nt,yzero((j-1)*3+1)-4,'UD','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
        end
        text(nt,yzero((j-1)*3+1)-0.2,num2str(abs(fit(1,1,rc(1)*(i-1)+j))),'fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
        text(nt,yzero((j-1)*3+1)-2.2,num2str(abs(fit(1,2,rc(1)*(i-1)+j))),'fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
        text(nt,yzero((j-1)*3+1)-4.2,num2str(abs(fit(1,3,rc(1)*(i-1)+j))),'fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top');
    end
    axis tight
    set(gca,'ycolor','w','ytick',[],'tickdir','out','xminortick','on',...
        'ylim',[-7,max(yzero)]+1,'ticklength',[0.02,0.02]);%,'grid','on')
    box off
    xlabel('Time(s)')
end
% return
% %     plot(t,ob0,'color',ones(1,3)*0.7,'linewidth',1.5);
% %     hold on
% %     plot(t,syn0,'k');
%     
%     %axis ij
%     if nt<0
%         plot([0,0],-[1,2*rc(1)+1],'k--');
%     end
%     set(gca,'xlim',[nt,mt],'ylim',[-2*rc(1)-1,-1]);
%     %axis off
%     
%     fit0=fit((i-1)*rc(1)+1:min(i*rc(1),end));
%     if ~isempty(station);station0=station((i-1)*rc(1)+1:min(i*rc(1),end),:);end
%     if ~isempty(phacom);phacom0=phacom((i-1)*rc(1)+1:min(i*rc(1),end),:);end
%     for j=1:size(ob0,2)
%         if ~isempty(station)
% %         text(nt+(mt-nt)*0.00,-2*j-0.025,[num2str(mm(j),2),'/',num2str(fit0(j),2)],...
% %             'fontsize',5,'HorizontalAlignment','left','VerticalAlignment','top');
%                 text(nt+(mt-nt)*0.00,-2*j-0.03,{num2str(mm(j),2);num2str(fit0(j),2)},...
%             'fontsize',6,'HorizontalAlignment','left','VerticalAlignment','top');
%         end
%         if ~isempty(phacom)
%         %text(mt,-2*j+0.2,[station0(j,:),'/',phacom0(j,:)],'fontsize',5,'HorizontalAlignment','right');
% %         text(nt+(mt-nt)*0.00,-2*j+0.025,[station0(j,:),'/',phacom0(j,:)],...
% %             'fontsize',5,'HorizontalAlignment','left','VerticalAlignment','bottom');
%                 text(nt+(mt-nt)*0.00,-2*j+0.03,{station0(j,:);phacom0(j,:)},...
%             'fontsize',6,'HorizontalAlignment','left','VerticalAlignment','bottom');
%         end
%     end
%     set(gca,'ycolor','w','ytick',[],'tickdir','out','xminortick','on');%,'grid','on')
%     box off
%     xlabel('Time(s)')
% end
%plot(ob(:,(i-1)*rc(1)+1:i*rc(1)),syn(:,(i-1)*rc(1)+1:i*rc(1)))

%==========================================================================
function [mout1,mout2,mm]=vecnor(min1,min2,rat)
%
if nargin==2;
    rat=1;
end
sm=size(min1);
mm1=max(abs(min1));mm2=max(abs(min2));
mm=max(mm1,mm2);
for i=1:sm(2)
    if mm(i)==0
        continue;
    end
    min1(:,i)=min1(:,i)/mm(i)*rat;%-i*2;
    min2(:,i)=min2(:,i)/mm(i)*rat;%-i*2;    
end
mout1=min1;mout2=min2;