function Ginsar=jt_sarG2(insar,locasub,dep,nsub,locasar,fault,gridsize,trup_con,sub_dex)
%


[slipg1]=slip_okada(locasub,dep,ones(nsub,1),locasar,fault,gridsize);


ss=size(slipg1);
slipsyn1=zeros(ss(1),ss(2));
for i=1:ss(2)
    slipsyn1(:,i)=slipg1(:,i,1).*insar(:,4)+slipg1(:,i,2).*insar(:,5)+slipg1(:,i,3).*insar(:,6);
end

Ginsar=geojt_conG(slipsyn1,trup_con,0);


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