function slip1=pltslip_pcolor(xx,yy,slip,maxslip,colorm,source,gridsize)
%


xx=mat2mat0(xx);
yy=mat2mat0(yy);


slip=mat2mat0(slip);
slip(2:end-1,[1,end])=slip(2:end-1,[2,end-1]);
slip([1,end],2:end-1)=slip([2,end-1],2:end-1);
slip(1,1)=(slip(1,2)+slip(2,1))/2;
slip(1,end)=(slip(1,end-1)+slip(2,end))/2;
slip(end,1)=(slip(end-1,1)+slip(end,2))/2;
slip(end,end)=(slip(end-1,end)+slip(end,end-1))/2;

slip1=slip;

bei=floor(sqrt(3e4/numel(slip)));% bei=3;

xx1=interm(xx,bei);
yy1=interm(yy,bei);

slip=interm(slip,bei);
slip(end,end)=maxslip;
%slip(slip<inter/2)=inter/2;
%pcolor(xx1,yy1,slip);
%shading interp;



pcolor_zh(xx1,yy1,slip);


%pcolor(xx1,yy1,slip);
axis ij
%imagesc(xx1(1,:),yy1(:,1),slip);
colormap(jet_zh(colorm,2));
%colorbar
%return

hold on

zbar=colorbar;
po_bar=get(zbar,'position');
set(get(zbar,'title'),'string','m');
% set(zbar,'position',[po_pic(1)+po_pic(3)./2,po_pic(2),po_pic(3)./2,po_pic(4)],'ticklength',[0.05,0.05]*5/5 ...
%     ,'fontname','times new roman','ylim',[0,maxslip],'ytick',0:maxslip/colorm:maxslip);

xyr=size(slip).*gridsize;
xyr=xyr(1)/xyr(2);
yshift=(po_bar(4)-po_bar(2))*xyr/1;

% set(zbar,'position',[po_bar(1)-po_bar(3)./2,po_bar(2)+yshift,po_bar(3)./2,po_bar(4)-2*yshift],'ticklength',[0.05,0.05]*5/5 ...
%     ,'fontname','times new roman','ylim',[0,maxslip],'ytick',0:maxslip/colorm:maxslip);
hold on
contour(xx1,yy1,slip,colorm-1,'color',ones(1,3)*1,'linewidth',1);%0.25

%contour(xx1,yy1,slip,colorm/10-1,'color',ones(1,3)*1,'linewidth',1.25);%0.25

%pltstar(8,[0,(source(1)-0.5)*gridsize(1)],gridsize(1)*3/3,0.48,'w','k',0.5);