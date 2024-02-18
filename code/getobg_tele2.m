function [ob,g,loca,mm,srate,leng,ga] = getobg_tele2(weight,ifrecal,obfile,epi,depth,fault,grid,gridsize,source,gpath,miu,band,isdis)
%(weight,ifrecal,obfile,depth,grid,gpath,miu,band,isdis)
% 获得远震记录的观测波形，计算格林函数，并进行滤波
if weight
    if ifrecal
        load(obfile);
        
%         nsta = size(ob,2);
%         dda=[repmat(depth,[size(loca,1),1]),da];
        [locasub,dep] = get_subloca(fault,grid([2,1]),gridsize([2,1]),[1,source([2,1])],loca,epi);
        dep=dep+depth-dep(source(1));
        nsta = size(loca,1);
        nsub = size(locasub,1);
        locasuba = repmat(locasub,[nsta,1]);
        depa = repmat(dep,[nsta,1]);
        locaa = repmat_zh(loca,nsub);
        daa = da_zh(locaa,locasuba,2);
        dda = [depa,daa];
        gt=zeros(size(ob,1),6,nsub*nsta);
        %change
        t0=round(get_time_new(da(:,1)/111.1949,depth,'P')*srate);
        t01=repmat(t0,[1,nsub]);
        t0=t01';
        t0=t0(:);%change
        tic
        for i=1:6
            M=zeros(6,1);M(i)=1;
            gtm = seekg_wang(gpath,dda,M,srate);
            %t=round(get_time_new(daa(:,1)/111.1949,depth,'P')*srate);
% 			t2=round(get_time_new(daa(:,1)/111.1949,dda(:,1),'P')*srate);
% 			t=diag(t2);
            for j=1:size(gtm,3)
                gt(:,i,j)=gtm(t0(j)-leng_before*srate:t0(j)+(leng-leng_before)*srate-1,1,j);
            end
            
        end
        gt=reshape(gt,size(ob,1),6,nsub,nsta);
        for i=1:nsub
            gt(:,:,i,:)=gt(:,:,i,:)*miu(i)/3e10;
        end
        g=gt;
        
%         if depth<100
%             t=round(get_time_new(da(:,1)/111.1949,depth,'P')*srate);
%         else
%             t=round(dbgrn_get_time(da(:,1)/111.1949,depth,'P')*srate);
%         end
%         g=zeros(650*srate,6,nsta);
%         for i=1:nsta
%             g0=gt(t(i)-leng_before*srate:t(i)+(650-leng_before)*srate-1,:,i);
% %             g0=gt(t(i)-10*srate:t(i)+640*srate-1,:,i);
%             g(1:size(g0,1),:,i)=g0;
%         end
%         g=g*miu/3e10;
%         g(size(ob,1)+1:end,:,:)=[];
        
%         gt = g;
%         for i = 1:prod(grid)-1
%             gt = cat(4,gt,g);
%         end
%         g = permute(gt,[1,2,4,3]);
        toc
        if isdis
            g = cumsum(g)/srate;
            ob = cumsum(ob)/srate;
        end
        mm = mm(:,1:4);
        disp('Wave and Green''s function of teleseismic record have been calculated!');
        save obg_tele.mat ob g loca mm srate leng
    else
        load('obg_tele.mat');
        disp('Wave and Green''s function of teleseismic record have been loaded!');
    end
    [bb,aa] = butter(4,band*2/srate);
    ob = filter(bb,aa,ob);
    g = filter(bb,aa,g);
else
    ob = [];g = [];loca = [];mm = [];srate = 0;leng = 0;
end
end
