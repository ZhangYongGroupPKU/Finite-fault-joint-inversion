function [ob,g,loca,mm,srate] = getobg_sm(weight,ifrecal,obfile,epi,depth,fault,grid,gridsize,source,gpath,miu,band,isdis)
% 获得强震记录的观测波形，计算格林函数，并进行滤波

if weight
    if ifrecal
        load(obfile);
        [locasub,dep] = get_subloca(fault,grid([2,1]),gridsize([2,1]),[1,source([2,1])],loca,epi);
        dep=dep+depth-dep(source(1));
        nsta = size(loca,1);
        nsub = size(locasub,1);
        locasuba = repmat(locasub,[nsta,1]);
        depa = repmat(dep,[nsta,1]);
        locaa = repmat_zh(loca,nsub);
        daa = da_zh(locaa,locasuba,2);
        dda = [depa,daa];
        tic
        bda = da_zh(locasuba,locaa,2);
        bfai = bda(:,2)+180;
        sfai = sin(bfai*pi/180);cfai = cos(bfai*pi/180);
        for i = 6:-1:1
            M = zeros(6,1);M(i) = 1;
            g0 = seekg_wang(gpath,dda,M,srate);
            g0(size(ob,1)+1:end,:,:) = [];
            for j = 1:size(g0,3)
                v = g0(:,1,j);r = g0(:,2,j);t = g0(:,3,j);
                g0(:,1,j) = r*sfai(j)-t*cfai(j); % E
                g0(:,2,j) = r*cfai(j)+t*sfai(j); % N
                g0(:,3,j) = v; % U
            end
            g1 = reshape(g0(:,1,:),[size(g0,1),nsub,nsta]);
            g2 = reshape(g0(:,2,:),[size(g0,1),nsub,nsta]);
            g3 = reshape(g0(:,3,:),[size(g0,1),nsub,nsta]);
            g(:,i,:,:) = cat(3,g1,g2,g3);
        end
        toc
        for i = 1:nsub
            g(:,:,i,:) = g(:,:,i,:)*miu(i)/3e10;
        end

        ob = permute(ob,[1,3,2]);
        ob = ob(:,:);
        if isdis == 0
        elseif isdis == 1
            g = cumsum(g)/srate;
            ob =cumsum(ob)/srate;
%             ob =cumsum(ob)/srate;%若输入为加速度，积分两次，如输入为速度则注释此行
        else
            error('Wrong ISDIS!!');
        end
        loca=repmat(loca,[3,1]);mm=repmat(mm,[3,1]);
        disp('Wave and Green''s function of strong-motion record have been calculated!');
        save obg_sm.mat ob g loca mm srate
    else
        load('obg_sm.mat');
        disp('Wave and Green''s function of strong-motion record have been loaded!');
    end
    [bb,aa] = butter(4,band*2/srate);
    ob = filter(bb,aa,ob);
    g = filter(bb,aa,g);
    %g = filter(bb,aa,g);
else
    ob = [];g = [];loca = [];mm = [];srate = 0;
end
