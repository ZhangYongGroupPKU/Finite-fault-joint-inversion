function [dex,dey]=substfs_plot(epi,grida,gridsize,source,substf,srate)
% This is a function to plot the sub-fault source time functions 
% 
% close all
substf=[substf;zeros(3,prod(grida))];
mastf=max(substf(:));
substf=substf/max(substf(:));
%substf=substf/1.916e18;
b=reshape(substf,size(substf,1),grida(1),grida(2));                       % substf�������������Ӷϲ�������������Ӷϲ����
if mod(grida(2),2)==0
    xz=(0:grida(2))-source(2)+1;xz=repmat(xz*gridsize(2),[grida(1),1]);xz=xz-gridsize(1)/2; % �Ѹ�20190823
else
    xz=(1:grida(2)+1)-source(2);xz=repmat(xz*gridsize(2),[grida(1),1])-gridsize(2)/2;
end
dex=[min(xz(:)),max(xz(:))];                                              % ������ķ�Χ
dey=[0,grida(1)*gridsize(1)];                                             % ������ķ�Χ

% ------ Plot the epi----
figure
plot(0,source(1)*gridsize(1)-gridsize(1)/2,'p','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','y','MarkerSize',20);% ����Դλ��
set(gca,'linewidth',2);                                                   % ��������Ӵ�
axis([dex,0,grida(1)*gridsize(1)]);                                       % x,y�ķ�Χ
set(gca,'xtick',dex(1):gridsize(1):dex(2),'FontSize',10)                  % ����x�������
set(gca,'ytick',dey(1):gridsize(1):dey(2),'FontSize',10);                 % ����y�������
set(gca,'YDir','reverse')                                                 % �����ᷭת�����µ������μ�С
set(gca,'dataaspectratio',[1,cosd(epi(1)),1])                             % �����ݺ���ı���
grid on
set(gca,'GridLineStyle','-');                                             % �������߱��ʵ��
xlabel('Distance along strike (km)','LineWidth',10,'FontSize',10);
ylabel('Distance down dip (km)','LineWidth',10,'FontSize',10);
hold on

% ----- Plot the source time functions ----
st=[];
dp=gridsize(1);                                                           % �������Ӷϲ�ĳ߶�
str=gridsize(2);                                                          % �������Ӷϲ�ĳ߶�
for j=1:grida(2)                                                          % ��ʾ�������Ӷϲ�ĸ���
    stf=b(:,:,j);
    for i=1:grida(1)                                                      % ��ʾ�������Ӷϲ�ĸ���
        y=i*dp-stf(:,i)*dp;
        nt=str/(size(substf,1)+srate*1);                                
        x=min(xz(:))+str*(j-1):nt:min(xz(:))+str+str*(j-1);               % x y ��γ�ȱ�����ͬ 23��Ҫ��
        x=x';
        x(size(y,1)+1:end,:)=[];
        sa=[x,y];sa=[sa;sa(1,:)];
        st=[st;sa];
        plot(x,y);
        fill(x,y,'c');
        hold on 
    end
end
title(['max amplitude ',num2str(mastf)]);
end
