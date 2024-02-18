function G=conG_wch(green,delt,fault,ndot,trup_con,inter,rake)
%       G=conG_wch(green,delt,fault,lensub,trup_con,inter,rake);
%==========================================================================
% G=conG_wch(green,delt,fault,grid,ndot,trup_con,inter,rake)
%--------------------------------------------------------------------------
%Input :green:[leng,6,nsub,nsta]
%
%Output:
%
%--------------------------------------------------------------------------
%                                 Zhang Yong,  2009/12/16 16:38
%               At Institute of Geophysics, China Earthquake Administration
%==========================================================================

% if nargin==7
%     rake=0;
% end

onet=find(trup_con~=0);
%zerot=find(trup_con==0);

sf=size(fault);
nsub=size(delt,1);
nsta=size(delt,2);
sg=size(green);%sg=[leng,6,nsub,nsta];
leng=sg(1);
if length(sg)==3
    % there is a path approximation from subfaults to stations.
else
    % path effect for each subfault to each station is considered!
end

if rake==1
    %sizemat=[grid(1,2),sum(grid(:,1))];
else
    g=green2g(green,fault(1,:));
    %sizemat=grid;
end

trup_con=reshape(trup_con,[nsub,ndot]);


if rake==1
%     leng
%     nsta
%     length(onet)
    G=zeros(leng*nsta,length(onet)*2);

    g1=zeros(sg(1),nsta,nsub);g2=zeros(sg(1),nsta,nsub);
    for ijk=1:nsub
        if length(size(green))==3
            g10=green2g(green,[fault(ijk,1:2),fault(ijk,3)-45]);
            g20=green2g(green,[fault(ijk,1:2),fault(ijk,3)+45]);
        elseif length(size(green))==4
            g10=green2g(green(:,:,ijk,:),[fault(ijk,1:2),fault(ijk,3)-45]);
            g20=green2g(green(:,:,ijk,:),[fault(ijk,1:2),fault(ijk,3)+45]);
        end
%         size(g10)
%         size(g1)
        g1(:,:,ijk)=g10(:,:);
        g2(:,:,ijk)=g20(:,:);
    end
    %numsub=[0;cumsum(grid(:,1))*grid(1,2)];

    for i=1:nsta
        kk=1;
        
        g11(:,:)=g1(:,i,:);
        g22(:,:)=g2(:,i,:);
%         subG1=con_subG(g11,delt(:,i),rake);
%         subG2=con_subG(g22,delt(:,i),rake);
        subG1=g11;
        subG2=g22;
        for j=1:ndot
%             zt=(trup_con(:,j)==0);
%             sG1=subG1;sG1(:,zt)=[];%subG1(:,zt)=[];
%             sG2=subG2;sG2(:,zt)=[];%subG2(:,zt)=[];

            ot=(trup_con(:,j)~=0);
            sG11=subG1(:,ot);
            sG22=subG2(:,ot);%每一stf时刻可以破裂的点的格林函数时间序列
            
            
            ss=size(sG11);%ss=size(subG1);
            
            %             G((i-1)*leng+(j-1)*inter+1:i*leng,kk:ss(2)+kk-1)=sG11(1:end-(j-1)*inter,:);
            %             G((i-1)*leng+(j-1)*inter+1:i*leng,(kk:ss(2)+kk-1)+length(onet))=sG22(1:end-(j-1)*inter,:);
            
            delt0=delt(ot,i);
            for ijk=1:numel(delt0)                
                tshift=(j-1)*inter+delt0(ijk); %old:tshift=(j-1)*inter;
                if tshift<0
                    % for teleseismic data inversion, it is always >0
                    % for near field data inversion, delt0==0
                    disp(num2str(tshift))
                    
                else
                    G((i-1)*leng+tshift+1:i*leng,kk+ijk-1)=sG11(1:end-tshift,ijk);
                    G((i-1)*leng+tshift+1:i*leng,kk+ijk-1+length(onet))=sG22(1:end-tshift,ijk);
                end
            end
            kk=kk+ss(2);
        end
    end
else
    G=zeros(leng*nsta,length(onet));%mem(G)
    for i=1:nsta
        kk=1;
        for j=1:ndot
            %subG=con_subG(g(:,:,i),delt(:,i));change
            subG=g(:,:,i);
            zt=(trup_con(:,j)==0);
            subG(:,zt)=[];

            ss=size(subG);

            G((i-1)*leng+(j-1)*inter+1:i*leng,kk:ss(2)+kk-1)=subG(1:end-(j-1)*inter,:);
            kk=kk+ss(2);
        end
    end
end
%==========================================================================

function subG=con_subG(g,delt,rake)
%
if nargin==2
    rake=0;
end

subG=zeros(size(g,1),length(delt));
if rake==1
    zvec=zeros(max(abs(delt(:))),1);
    for i=1:length(delt)
        if delt(i)>0
            %zg=[zeros(delt(i),1);g(1:end-delt(i))];
            zg=[zvec(1:delt(i));g(1:end-delt(i),i)];
        elseif delt(i)<0
            %zg=[g(1-delt(i):end);zeros(-delt(i),1)];
            zg=[g(1-delt(i):end,i);zvec(1:-delt(i))];
        else
            zg=g(:,i);
        end
        subG(:,i)=zg;
    end
else% rake==0
    zvec=zeros(max(abs(delt(:))),1);
    for i=1:length(delt)
        if delt(i)>0
            %zg=[zeros(delt(i),1);g(1:end-delt(i))];
            zg=[zvec(1:delt(i));g(1:end-delt(i))];
        elseif delt(i)<0
            %zg=[g(1-delt(i):end);zeros(-delt(i),1)];
            zg=[g(1-delt(i):end);zvec(1:-delt(i))];
        else
            zg=g;
        end
        subG(:,i)=zg;
    end
end

return

%===================================end====================================
% if nargin==2
%     subG=zeros(length(g),length(delt));
%     for i=1:length(delt)
%         if delt(i)>0
%             zg=[zeros(delt(i),1);g(1:end-delt(i))];
%         elseif delt(i)<0
%             zg=[g(1-delt(i):end);zeros(-delt(i),1)];
%         else
%             zg=g;
%         end
%         subG(:,i)=zg;
%     end
% else
%     % rake==1
%     subG=zeros(length(g),length(delt));
%     for j=1:length(delt)
%         if delt(j)>0
%             zg=[zeros(delt(j),1);g(1:end-delt(j),j)];
%         elseif delt(j)<0
%             zg=[g(1-delt(j):end,j);zeros(-delt(j),1)];
%         else
%             zg=g(:,j);
%         end
%         subG(:,j)=zg;
%     end
% end