function [dex,dey]=substfs_plot(epi,grida,gridsize,source,substf,srate)
% This is a function to plot the sub-fault source time functions 
% 
% close all
substf=[substf;zeros(3,prod(grida))];
mastf=max(substf(:));
substf=substf/max(substf(:));
%substf=substf/1.916e18;
b=reshape(substf,size(substf,1),grida(1),grida(2));                       % substf行数×倾向方向子断层个数×走向方向子断层个数
if mod(grida(2),2)==0
    xz=(0:grida(2))-source(2)+1;xz=repmat(xz*gridsize(2),[grida(1),1]);xz=xz-gridsize(1)/2; % 已改20190823
else
    xz=(1:grida(2)+1)-source(2);xz=repmat(xz*gridsize(2),[grida(1),1])-gridsize(2)/2;
end
dex=[min(xz(:)),max(xz(:))];                                              % 走向方向的范围
dey=[0,grida(1)*gridsize(1)];                                             % 倾向方向的范围

% ------ Plot the epi----
figure
plot(0,source(1)*gridsize(1)-gridsize(1)/2,'p','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','y','MarkerSize',20);% 画震源位置
set(gca,'linewidth',2);                                                   % 将坐标轴加粗
axis([dex,0,grida(1)*gridsize(1)]);                                       % x,y的范围
set(gca,'xtick',dex(1):gridsize(1):dex(2),'FontSize',10)                  % 设置x轴的坐标
set(gca,'ytick',dey(1):gridsize(1):dey(2),'FontSize',10);                 % 设置y轴的坐标
set(gca,'YDir','reverse')                                                 % 将纵轴翻转，从下到上依次减小
set(gca,'dataaspectratio',[1,cosd(epi(1)),1])                             % 设置纵横轴的比例
grid on
set(gca,'GridLineStyle','-');                                             % 将网格线变成实线
xlabel('Distance along strike (km)','LineWidth',10,'FontSize',10);
ylabel('Distance down dip (km)','LineWidth',10,'FontSize',10);
hold on

% ----- Plot the source time functions ----
st=[];
dp=gridsize(1);                                                           % 倾向方向子断层的尺度
str=gridsize(2);                                                          % 走向方向子断层的尺度
for j=1:grida(2)                                                          % 表示走向方向子断层的个数
    stf=b(:,:,j);
    for i=1:grida(1)                                                      % 表示倾向方向子断层的个数
        y=i*dp-stf(:,i)*dp;
        nt=str/(size(substf,1)+srate*1);                                
        x=min(xz(:))+str*(j-1):nt:min(xz(:))+str+str*(j-1);               % x y 的纬度必须相同 23需要改
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
