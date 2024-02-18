function [Ggps,G_gps]=jt_gpsG(locagps,locasub,locadep,fault,gridsize,trup_con)
%==========================================================================
%       Ginsar=jt_gpsG(locagps,locasub,locadep,fault,gridsize,trup_con)
%    Here is a function to construct the matrix for GPS data in joint
%    inversion
%--------------------------------------------------------------------------
%    Input
%           locagps: location of observed GPS stations.[lat,long]
%           locasub: location of subfaults
%           locadep: depth of subfault
%             fault: fault parameters of each subfault
%          gridsize: size of subfault
%          trup_con: smpling dots allowed to slip in joint inversion under
%          the constraints of rupture speed and source duration 
%    Output
%           G_insar: matrix result. with the size of
%                    [nsub,length(trup_con)]
%--------------------------------------------------------------------------
%                               Zhang Yong,  2009/12/22 09:51
%                 Institute of Geophysics, China Earthquake Administration
%==========================================================================

fault(:,3)=fault(:,3)-45;
[slipg1]=slip_okada(locasub,locadep,ones(length(locadep),1),locagps,fault,gridsize);
fault(:,3)=fault(:,3)+90;
[slipg2]=slip_okada(locasub,locadep,ones(length(locadep),1),locagps,fault,gridsize);


% fault(:,3)=fault(:,3)-45;
% [slipg1]=slip_okada_parea(locasub,locadep,ones(length(locadep),1),locagps,fault,gridsize);
% fault(:,3)=fault(:,3)+90;
% [slipg2]=slip_okada_parea(locasub,locadep,ones(length(locadep),1),locagps,fault,gridsize);


%slipg1 and slipg2 have the size of [nsta,nsub,3]

% ss=size(slipg1);
% slipsyn1=zeros(ss(1),ss(2));
% for i=1:ss(2)
%     slipsyn1(:,i)=slipg1(:,i,1).*insar(:,4)+slipg1(:,i,2).*insar(:,5)+slipg1(:,i,3).*insar(:,6);
% end

ss=size(slipg1);
G_gps=zeros(ss(1)*3,ss(2)*2);
for i=1:3
    G_gps((1:ss(1))+(i-1)*ss(1),1:ss(2))=slipg1(:,:,i);
    G_gps((1:ss(1))+(i-1)*ss(1),ss(2)+(1:ss(2)))=slipg2(:,:,i);
end

rake=1;% since there are slipg1 and slipg2, so the rake variations are considered
Ggps=geojt_conG(G_gps,trup_con,rake);
%Ginsar=geojt_conG([slipsyn1,slipsyn2],trup_con,1);

%=================================end======================================



% 
% 
% obs=gps(:,3:5);
% 
% sizemat=[grid(1,2),sum(grid(:,1))];sizemat1=[grid1(1,2),sum(grid1(:,1))];
% smG=sm_space_ad(sizemat,1);smG1=sm_space_ad(sizemat1,1);
% ss=size(smG);ss1=size(smG1);
% sparG=[smG,sparse(ss(1),ss1(2));sparse(ss1(2),ss(1)),smG1];
% 
% ss=size(sparG);
% sparG=[sparG,sparse(ss(1),ss(2));sparse(ss(1),ss(2)),sparG];
% %[X,rsq]=cgls_zh0(G_gps,obs(:),500,'X(X<0)=0;');
% 
% lambda=7e-2;
% [X,rsq]=cgls_zh0_spar(G_gps,lambda.*sparG,[obs(:);zeros(840,1)],500,'X(X<0)=0;');