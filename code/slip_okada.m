function [slipend]=slip_okada(locasub,dep,slip,loca,fault,gridsize)
%===============================================================================
%  slipend=slip_chen(sizemat,locasub,slip,loca,fault,grid);
%  Here it is a function for static slip calculation on ground surface
%-------------------------------------------------------------------------------
%  Input:
%   locasub:  location of subfaults.structure as:[latitude(:),longitude(:)] 
%       dep:  depth of all subfaults (Units:km)
%      slip:  slip amplitudes of all subfaults
%      loca:  location of stations or observation points
%     fault:  fault parameters of the fault plane. [strike,dip,rake]. while
%             its size(fault,1) more than one, each row corresponds a
%             mechanism of each subfault, then it has a size of [nsta,3] 
%  gridsize:  'grid' describle the size of this
%             subfault as [size in dip direction, size in strike direction] with
%             the units of km
%
%  Output:
%   slipend:  if slip(1)==1, it has the size:[nsta,nsub,3], correspond components
%              E-W,N-S,U-D. elseif slip(1)~=1, it has size of [nsta,3]
%         
%-------------------------------------------------------------------------------
%   For more details, please see Okada (1985) about the theoretical back 
%   ground.
%
%                                   Zhang Yong,  Chen Yun-Tai,  Xu Li-Sheng
%                                   2009/06/02  18:43     CEAIGP
%===============================================================================

sl=size(locasub);
nsub=sl(1);
sl=size(loca);
nsta=sl(1);

sf=size(fault);
if sf(1)==1
    fault=repmat(fault,[nsub,1]);
end
sg=size(gridsize);
if sg(1)==1
    gridsize=repmat(gridsize,[nsub,1]);
end

slipend=zeros(nsta,nsub,3);
for i=1:nsub
    if slip(i)==0
        continue;
    end
    strike=fault(i,1);
    dip=fault(i,2);
    rake=fault(i,3);
    slipvec=slip(i);
    DEP=dep(i);
    LL=gridsize(i,2);
    WW=gridsize(i,1);
    
    
    loca0=locasub(i,:);
    [d,a]=distance(loca0,loca);
    d=d*6371*pi/180;
    %[da]=distazim_zh(loca,loca0);d=da(:,1);a=da(:,2);
        
    
    %X=d.*sind(180-(a-strike+90));Y=d.*cosd(180-(a-strike+90));
    X=d.*sind(a-strike+90);Y=d.*cosd(a-strike+90);
    
    
    DISL1=slipvec*cosd(rake);
    DISL2=slipvec*sind(rake);
    DISL3=0;
    
%     strike=zmat(i,6);
%     DISL11=DISL1*sind(strike)+DISL2*cosd(dip)*sind(strike);
%     DISL11=DISL1*sind(strike)+DISL2*cosd(dip)*sind(strike);

    ALP=0.5;
    AL1=-LL/2;
    AL2=LL/2;
    AW1=-WW/2;
    AW2=WW/2;
    
    SD=sind(dip);
    CD=cosd(dip);


    OU=okada2d(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,...
        SD,CD,DISL1,DISL2,DISL3);
    
    slipend(:,i,1)=OU{1};
    slipend(:,i,2)=OU{2};
    slipend(:,i,3)=OU{3};
end

% strike=repmat(fault(:,1)',[nsta,1,1]);
% slipx=slipend(:,:,1).*sind(strike)+slipend(:,:,2).*sind(strike-90);
% slipy=slipend(:,:,1).*cosd(strike)+slipend(:,:,2).*cosd(strike-90);
% slipend(:,:,1)=slipx;
% slipend(:,:,2)=slipy;

strike=fault(:,1);
for i=1:nsub
    slipx=slipend(:,i,1)*sind(strike(i))+slipend(:,i,2)*sind(strike(i)-90);
    slipy=slipend(:,i,1)*cosd(strike(i))+slipend(:,i,2)*cosd(strike(i)-90);
    slipend(:,i,1)=slipx;
    slipend(:,i,2)=slipy;
end

if slip(1)~=1
    slipend1(:,:)=sum(slipend,2);
    slipend=slipend1;
end
return
%=====================================end=======================================