function rate=rup_pltslipshot2(X,sizemat,lensub,para,nap,maxslip,numth,srate)
%

if sizemat(1)==1
    X=reshape(X,[lensub,sizemat(1)*sizemat(2)]);
    X=[X;X;X];
    sizemat(1)=3;
    %     sizegrid(1)=sizegrid(1)./3;
    %     epi(1)=2;
elseif sizemat(2)==1
    X=reshape(X,[lensub,sizemat(1)*sizemat(2)]);
    X=[X,X,X];
    sizemat(2)=3;
    %     sizegrid(2)=sizegrid(2)./3;
    %     epi(2)=2;
end

if nargin<6
    if nap==1
        X1=reshape(X,[lensub,sizemat(1),sizemat(2)]);
        X2=zeros([ceil(lensub./para(1)),sizemat(1),sizemat(2)]);
        for j=1:ceil(lensub./para(1))
            if para(1)==1
                X2(j,:,:)=(X1((j-1)*para(1)+1:min(j*para(1),end),:,:));
            else
                X2(j,:,:)=sum(X1((j-1)*para(1)+1:min(j*para(1),end),:,:));
            end
        end
        maxslip=max(X2(:));
    else
        maxslip=max(sum(reshape(X(:),[lensub,length(X(:))./lensub])));
    end
end

len=prod(para);
hang=para(2);
lie=para(3);
X=reshape(X,[lensub,sizemat(1),sizemat(2)]);

sor=zeros(hang,lie);
for i=1:hang
    for j=1:lie
        sor(i,j)=(i-1)*lie+j;
    end
end
figure
xz=zeros(sizemat);
nx=0;

rate=[];
for i=para(1):para(1):min(len,lensub)
    nx=nx+1;
    subplot(hang,lie,sor(nx))
    if nap==1
        xz(:,:)=sum(X(i-para(1)+1:i,:,:),1);
        nz=flipud(interm(xz,30));
        if max(nz(:))==maxslip
            rate=nz;
        end
        nz(end,end)=maxslip;

        %nz(end,end)=0.7;
    else
        xz(:,:)=sum(X(1:i,:,:),1);

        slip=xz;
        slip=mat2mat0(slip);
        slip(2:end-1,[1,end])=slip(2:end-1,[2,end-1]);
        slip([1,end],2:end-1)=slip([2,end-1],2:end-1);
        slip(1,1)=(slip(1,2)+slip(2,1))/2;
        slip(1,end)=(slip(1,end-1)+slip(2,end))/2;
        slip(end,1)=(slip(end-1,1)+slip(end,2))/2;
        slip(end,end)=(slip(end-1,end)+slip(end,end-1))/2;
        xz1=slip;


        nz=flipud(interm(xz1,6));
        nz(end,end)=maxslip;

    end
    %set(gca,'layer','top','linewidth',1);


    pcolor(nz);
    hold on
    shading interp
    %colormap(colorprod12([0 0.28 0.30 0.5 0.7 1.2 1.45 1.5 1.75 1.9 2.95 2.95]));
    %colormap(colorprod12([0 0.36 0.63 0.85 1.0 1.2 1.5 1.82 2 2.2 2.65 2.95]));
    %colormap(colorprod12([0 0.57 0.76 0.95 1.1 1.4 1.85 2 2.1 2.2 2.65 2.95]));
    %     colormap(colorprod12([0 0.39 0.76 0.95 1.1 1.4 1.85 2 2.1 2.2 2.65 2.95]));
    %     colormap(colorprod12([0 0.1 0.35 0.5 0.7 1.1 1.4 1.6 1.75 1.9 2.65 2.95]));

    %rate
    %colormap(colorprod12([0 0.29 0.45 0.7 1.1 1.5 1.8 2.0 2.1 2.2 2.45 2.95]));
    %colormap(colorprod12([0 0.001 0.33 0.6 1.0 1.3 1.6 1.82 1.95 2.02 2.5 2.95]));

    %    axis equal

if nargin<7
numth=10;
else
end
    %num=7;
    colormap(jet_zh(numth,2));

    axis tight
    box on
    hold on
    contour(nz,numth-1,'color',[1,1,1]);
    set(gca,'xtick',[],'ytick',[],'linewidth',1,'layer','top');
    %ceil(max(slip(:)))-1
    %     if i==para(1)
    %         snz=size(nz);
    %         plot(snz(2)*2/3,snz(1)*2/3,'kh','markerfacecolor','w','markersize',18);
    %         %plot(snz(2)*4/15,snz(1)*5/6,'kh','markerfacecolor','w','markersize',18);%8.4
    %         %plot(snz(2)*11/35,snz(1)*9/12,'kh','markerfacecolor','w','markersize',18);%7.9
    %     end
end