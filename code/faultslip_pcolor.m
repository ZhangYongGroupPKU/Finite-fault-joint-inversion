function [slip1,start,stop]=faultslip_pcolor(grid,gridsize,source,slip,maxslip,numth,inter,alen)

start=[];
stop=[];

ss=size(slip);
if length(ss)>2
    slipstrike=slip(:,:,1);
    slipdip=slip(:,:,2);
    slip=sqrt(slipstrike.^2+slipdip.^2);
end

%xz=[1:grid(2)]-source(2);xz=repmat(xz*gridsize(2)-gridsize(2)/2,[grid(1),1]);
xz=[1:grid(2)]-source(2);xz=repmat(xz*gridsize(2),[grid(1),1]);
yz=[1:grid(1)]';yz=repmat(yz*gridsize(1)-gridsize(1)/2,[1,grid(2)]);

%pcolor(xz,yz,slip)
slip1=pltslip_pcolor(xz,yz,slip,maxslip,numth,source,gridsize);
hold on

if length(ss)>2
    dex=slip<inter;
    if nargin<8
    [start,stop]=arrow_zh_new(xz,yz,slipstrike,slipdip,[3,4]*(min(gridsize)/10).^0.5,dex);
    else
    [start,stop]=arrow_zh_new(xz,yz,slipstrike,slipdip,alen,dex);
    end
end

%aniseStar( [0,(source(1)-0.5)*gridsize(1)], gridsize(1)/2, 1, 'w');
pltstar(8,[0,(source(1)-0.5)*gridsize(1)],gridsize(1)/3,0.48,'w','k',1);

xlabel('Distance along strike (km)','fontname','times new roman','fontsize',16);
ylabel('Distance down dip (km)','fontname','times new roman','fontsize',16);

% set(gca,'linewidth',0.5,'ticklength',[0.02,0.02],'tickdir','out','layer','top','fontname','times new roman');
%set(gca,'linewidth',0.7,'ticklength',[0.02,0.02],'tickdir','out','ylim',[0,240],'layer','bottom','fontname','times new roman');
%==========================================================================
%==========================================================================
%==========================================================================
function [start,stop]=arrow_zh_new(x,y,slipx,slipy,len,dex)
%
[u,v]=quiverc_repair(x,y,slipx,slipy,len(1));

start=[u(1,:);v(1,:)];start=start';
stop=[u(2,:);v(2,:)];stop=stop';

start(dex,:)=[];stop(dex,:)=[];


for i=1:size(start,1);
[h,dot]=arrow0(start(i,:),stop(i,:),0.2,0.01,[20,50],ones(1,3)*0.5,ones(1,3)*0.5);
end

%arrow('start',start','stop',stop','tipangle',15,'baseangle',50,'length',len(2),'width',0.03,'facecolor',ones(1,3).*0.4,'edgecolor',ones(1,3).*0.4);


 function [uu,vv] = quiverc_repair(varargin)
% Modified version of Quiver to plots
% velocity vectors as arrows 
% with components (u,v) at the points (x,y) using the current colormap 

% Bertrand Dano 3-3-03
% Copyright 1984-2002 The MathWorks, Inc. 

%QUIVERC Quiver color plot.
%   QUIVERC(X,Y,U,V) plots velocity vectors as arrows with components (u,v)
%   at the points (x,y).  The matrices X,Y,U,V must all be the same size
%   and contain corresponding position and velocity components (X and Y
%   can also be vectors to specify a uniform grid).  QUIVER automatically
%   scales the arrows to fit within the grid.
%
%   QUIVERC(U,V) plots velocity vectors at equally spaced points in
%   the x-y plane.
%
%   QUIVERC(U,V,S) or QUIVER(X,Y,U,V,S) automatically scales the 
%   arrows to fit within the grid and then stretches them by S.  Use
%   S=0 to plot the arrows without the automatic scaling.
%
%   QUIVERC(...,LINESPEC) uses the plot linestyle specified for
%   the velocity vectors.  Any marker in LINESPEC is drawn at the base
%   instead of an arrow on the tip.  Use a marker of '.' to specify
%   no marker at all.  See PLOT for other possibilities.
%
%   QUIVERC(...,'filled') fills any markers specified.
%
%   H = QUIVERC(...) returns a vector of line handles.
%
%   Example:
%      [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%      z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
%      contour(x,y,z), hold on
%      quiverc(x,y,px,py), hold off, axis image
%
%   See also FEATHER, QUIVER3, PLOT. 
%   Clay M. Thompson 3-3-94
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.21 $  $Date: 2002/06/05 20:05:16 $ 
%-------------------------------------------------------------

%set(gca, 'color', 'blue');
% Arrow head parameters
alpha = 0.33; % Size of arrow head relative to the length of the vector
beta = 0.23;  % Width of the base of the arrow head relative to the length
autoscale = 0; % Autoscale if ~= 0 then scale by this.
plotarrows = 0; % Plot arrows
sym = '';

filled = 0;
ls = '-';
ms = '';
col = '';
lw=1;

nin = nargin;
% Parse the string inputs
while isstr(varargin{nin}),
  vv = varargin{nin};
  if ~isempty(vv) & strcmp(lower(vv(1)),'f')
    filled = 1;
    nin = nin-1;
  else
    [l,c,m,msg] = colstyle(vv);
    if ~isempty(msg), 
      error(sprintf('Unknown option "%s".',vv));
    end
    if ~isempty(l), ls = l; end
    if ~isempty(c), col = c; end
    if ~isempty(m), ms = m; plotarrows = 0; end
    if isequal(m,'.'), ms = ''; end % Don't plot '.'
    nin = nin-1;
  end
end


error(nargchk(2,5,nin));

% Check numeric input arguments
if nin<4, % quiver(u,v) or quiver(u,v,s)
  [msg,x,y,u,v] = xyzchk(varargin{1:2});
else
  [msg,x,y,u,v] = xyzchk(varargin{1:4});
end
if ~isempty(msg), error(msg); end

if nin==3 | nin==5, % quiver(u,v,s) or quiver(x,y,u,v,s)
  autoscale = varargin{nin};
end

% Scalar expand u,v
if prod(size(u))==1, u = u(ones(size(x))); end
if prod(size(v))==1, v = v(ones(size(u))); end

if autoscale,
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
  delx = diff([min(x(:)) max(x(:))])/n;
  dely = diff([min(y(:)) max(y(:))])/m;
  len = sqrt((u.^2 + v.^2)/(delx.^2 + dely.^2));
  autoscale = autoscale*0.9 / max(len(:));
  u = u*autoscale; v = v*autoscale;
end

%----------------------------------------------
% Define colormap 
% vr=sqrt(u.^2+v.^2);
% vrn=round(vr/max(vr(:))*64);
% CC=colormap;
% ax = newplot;
% next = lower(get(ax,'NextPlot'));
% hold_state = ishold;

%----------------------------------------------
% Make velocity vectors and plot them

x = x(:).';y = y(:).';
u = u(:).';v = v(:).';
%vrn=vrn(:).';
uu = [x;x+u];
vv = [y;y+v];
% uu = [x;x+u;repmat(NaN,size(u))];
% vv = [y;y+v;repmat(NaN,size(u))];
% vrn1= [vrn;repmat(NaN,size(u));repmat(NaN,size(u))];
% 
% uui=uu(:);  vvi=vv(:);  vrn1=vrn1(:); imax=size(uui);
% hold on
% 
%  for i=  1:3:imax-1
%     ii=int8(round(vrn1(i)));
%     if ii==0; ii=1; end        
%     c1= CC(ii,1);    c2= CC(ii,2);    c3= CC(ii,3);
%     plot(uui(i:i+1),vvi(i:i+1),'linewidth',lw,'color',[c1 c2 c3]);
% end
% 
% %----------------------------------------------
% % Make arrow heads and plot them
% if plotarrows,
%  
%   hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
%         x+u-alpha*(u-beta*(v+eps));repmat(NaN,size(u))];
%   hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
%         y+v-alpha*(v+beta*(u+eps));repmat(NaN,size(v))];
%   vrn2= [vrn;vrn;vrn;vrn];
% 
%  uui=hu(:);  vvi=hv(:);  vrn2=vrn2(:); imax=size(uui);
% 
%  for i=  1:imax-1
%     ii=int8(round(vrn2(i)));
%     if ii==0; ii=1; end   
%     c1= CC(ii,1);    c2= CC(ii,2);    c3= CC(ii,3);
%     plot(uui(i:i+1),vvi(i:i+1),'linewidth',lw,'color',[c1 c2 c3]);
%  end
% 
% else
%   h2 = [];
% end
% %----------------------------------------------
% 
% if ~isempty(ms), % Plot marker on base
%   hu = x; hv = y;
%   hold on
%   h3 = plot(hu(:),hv(:),[col ms]);
%   if filled, set(h3,'markerfacecolor',get(h1,'color')); end
% else
%   h3 = [];
% end
% 
% if ~hold_state, hold off, view(2); set(ax,'NextPlot',next); end
% 
% % if nargout>0, hh = [h1;h2;h3]; end
% % set(gca, 'color', [0 0 0],'Xcolor','w','Ycolor','w');
% % set(gcf, 'color', [0 0 0]);
% %set(gcf, 'InvertHardCopy', 'off');
% 
