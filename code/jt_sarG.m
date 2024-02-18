function Ginsar=jt_sarG(insar,locasub,dep,nsub,locasar,fault,gridsize,trup_con,sub_dex)
%


fault(:,3)=fault(:,3)-45;
[slipg1]=slip_okada(locasub,dep,ones(nsub,1),locasar,fault,gridsize);
fault(:,3)=fault(:,3)+90;
[slipg2]=slip_okada(locasub,dep,ones(nsub,1),locasar,fault,gridsize);

ss=size(slipg1);
slipsyn1=zeros(ss(1),ss(2));
for i=1:ss(2)
    slipsyn1(:,i)=slipg1(:,i,1).*insar(:,4)+slipg1(:,i,2).*insar(:,5)+slipg1(:,i,3).*insar(:,6);
end

slipsyn2=zeros(ss(1),ss(2));
for i=1:ss(2)
    slipsyn2(:,i)=slipg2(:,i,1).*insar(:,4)+slipg2(:,i,2).*insar(:,5)+slipg2(:,i,3).*insar(:,6);
end

Ginsar=geojt_conG([slipsyn1,slipsyn2],trup_con,1);


%--------------------------------------
if nargin>8
    sg=size(Ginsar);
    one_sub=zeros(sg(1),length(sub_dex));
    sumdex=[0;cumsum(sub_dex(:))];
    for i=1:length(sub_dex)
        one_sub(sumdex(i)+1:sumdex(i+1),i)=1;
    end
    
    Ginsar=[Ginsar,one_sub];
end
%=================================end======================================