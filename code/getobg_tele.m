function [ob,g,loca,mm,srate,leng,ga] = getobg_tele(weight,ifrecal,obfile,depth,grid,gpath,miu,band,isdis)
% 获得远震记录的观测波形，计算格林函数，并进行滤波
if weight
    if ifrecal
        load(obfile);
        nsta = size(ob,2);
        dda=[repmat(depth,[size(loca,1),1]),da];
        gt=zeros(650*srate,6,nsta);
        tic
        for i=1:6
            M=zeros(6,1);M(i)=1;
            [g0,depth]=seekg_wang(gpath,dda(:,:),M,srate);
            gt(1:size(g0,1),i,:)=g0(:,1,:);
        end
        depth=depth(1);
        
        ga = gt;
        
        if depth<100
            t=round(get_time_new(da(:,1)/111.1949,depth,'P')*srate);
        else
            t=round(dbgrn_get_time(da(:,1)/111.1949,depth,'P')*srate);
        end
        g=zeros(650*srate,6,nsta);
        for i=1:nsta
            g0=gt(t(i)-leng_before*srate:t(i)+(650-leng_before)*srate-1,:,i);
%             g0=gt(t(i)-10*srate:t(i)+640*srate-1,:,i);
            g(1:size(g0,1),:,i)=g0;
        end
        g=g*miu/3e10;
        g(size(ob,1)+1:end,:,:)=[];
        
        gt = g;
        for i = 1:prod(grid)-1
            gt = cat(4,gt,g);
        end
        g = permute(gt,[1,2,4,3]);
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
    [bb,aa] = butter(3,band*2/srate);
    ob = filter(bb,aa,ob);
    g = filter(bb,aa,g);
else
    ob = [];g = [];loca = [];mm = [];srate = 0;leng = 0;
end
end
