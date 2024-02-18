function pltobsyn_new(ob,syn,t,rc,rat,station,phacom,flag)
%==========================================================================

%==========================================================================

if nargin==4
    rat=0.05;
    station=[];
    phacom=[];flag=1;
elseif nargin==5
    station=[];
    phacom=[];flag=1;
elseif nargin==6
    phacom=[];flag=1;
elseif nargin==7
    flag=1;
end

wid0=1/(rc(2)*(1+rat)+rat);
wid=wid0*rat;
%high=1/rc(1);


fit=gfit0(ob,syn);
figure('color','w')

nt=min(t);mt=max(t);
for i=1:rc(2)
    axes('position',[(i-1)*wid0+wid*i,0.1,wid0,0.9]);
    
    dex0=(i-1)*rc(1)+1;
    if dex0>size(ob,2);
        axis off
        return;
    end
    
    ob0=ob(:,(i-1)*rc(1)+1:min(i*rc(1),end));
    syn0=syn(:,(i-1)*rc(1)+1:min(i*rc(1),end));
    if flag==1
        [ob0,syn0,mm]=vecnor(ob0,syn0,1-rat);
    else
        nsta0=size(ob0,2);
        mmob=1./max(abs(ob0));
        ob0=ob0*sparse(1:nsta0,1:nsta0,mmob)*(1-2*rat)+0.2;
        ob0=ob0+repmat((-2:-2:-2*nsta0),size(ob0,1),1);
        mmsyn=1./max(abs(syn0));
        syn0=syn0*sparse(1:nsta0,1:nsta0,mmsyn)*(1-2*rat)-0.2;
        syn0=syn0+repmat((-2:-2:-2*nsta0),size(ob0,1),1);
        mm=mmob;
    end
    
    plot(t,ob0,'k','linewidth',1.5);
    hold on
    plot(t,syn0,'r','linewidth',1.5);
    
    
%     plot(t,ob0,'color',ones(1,3)*0.7,'linewidth',1.5);
%     hold on
%     plot(t,syn0,'k');
    
    %axis ij
    % plot the line for x=0
    if nt<0
        %plot([0,0],-[1,2*rc(1)+1],'k--');
    end
    
    set(gca,'xlim',[nt,mt],'ylim',[-2*rc(1)-1,-1]);
    %axis off
    
    fit0=fit((i-1)*rc(1)+1:min(i*rc(1),end));
    
    if ~isempty(station);station0=station((i-1)*rc(1)+1:min(i*rc(1),end),:);end
    if ~isempty(phacom);phacom0=phacom((i-1)*rc(1)+1:min(i*rc(1),end),:);end
    for j=1:size(ob0,2)
        if ~isempty(station)
%         text(nt+(mt-nt)*0.00,-2*j-0.025,[num2str(mm(j),2),'/',num2str(fit0(j),2)],...
%             'fontsize',5,'HorizontalAlignment','left','VerticalAlignment','top');
                text(nt+(mt-nt)*0.00,-2*j-0.03,{num2str_zh(mm(j),2);num2str_zh(fit0(j),2)},...
            'fontsize',9,'HorizontalAlignment','left','VerticalAlignment','top','fontname','times');
        end
        if ~isempty(phacom)
        %text(mt,-2*j+0.2,[station0(j,:),'/',phacom0(j,:)],'fontsize',5,'HorizontalAlignment','right');
%         text(nt+(mt-nt)*0.00,-2*j+0.025,[station0(j,:),'/',phacom0(j,:)],...
%             'fontsize',5,'HorizontalAlignment','left','VerticalAlignment','bottom');
                text(nt+(mt-nt)*0.00,-2*j+0.03,{station0(j,:);phacom0(j,:)},...
            'fontsize',9,'HorizontalAlignment','left','VerticalAlignment','bottom','fontname','times');
        end
    end
    set(gca,'ycolor','w','ytick',[],'tickdir','out','xminortick','on','fontname','times');%,'grid','on')
    box off
    xlabel('Time(s)')
end
%plot(ob(:,(i-1)*rc(1)+1:i*rc(1)),syn(:,(i-1)*rc(1)+1:i*rc(1)))

%==========================================================================
function [mout1,mout2,mm]=vecnor(min1,min2,rat)
%
sm=size(min1);
mm1=max(abs(min1));mm2=max(abs(min2));
mm=max(mm1,mm2);
for i=1:sm(2)
    if mm(i)==0
        continue;
    end
    min1(:,i)=min1(:,i)/mm(i)*rat-i*2;
    min2(:,i)=min2(:,i)/mm(i)*rat-i*2;    
end
mout1=min1;mout2=min2;