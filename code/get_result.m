function [slip,stf,slipx,slipy,slipe,slipn,zslip]=get_result(X,lensub,sizemat,gridsize,fault,miu,srate)
%
if nargin<6
    miu=3e10;
end

lensub = lensub*srate;

nsub=sum(prod(sizemat'));
if numel(X)>nsub*lensub
    X1=sqrt(X(:,1).^2+X(:,2).^2);

    stf_sub=reshape(X1,[lensub,nsub]);
    %stf_sub=stf_sub*diag(miu)/3e10;%/prod(gridsize)/1e6;
    stf=sum(stf_sub/3e10*diag(miu(:)),2)*srate;
    %stf=sum(sum(reshape(X1,[lensub,sizemat]),2),3);
    
    Xslip=X/3e16/prod(gridsize);
    slip1(:,:)=sum(reshape(Xslip(:,1),[lensub,sizemat]));
    slip2(:,:)=sum(reshape(Xslip(:,2),[lensub,sizemat]));
    slip=sqrt(slip1.^2+slip2.^2);
    
%     slip=slipnf(sizemat,slip);
    
   
    if nargin>4
        slipx=slip1.*cosd(fault(:,3)-45)+slip2.*cosd(fault(:,3)+45);
        slipy=slip1.*sind(fault(:,3)-45)+slip2.*sind(fault(:,3)+45);
        
        slipe=slipx.*sind(fault(:,1))+slipy.*sind(fault(:,1)-90);
        slipn=slipx.*cosd(fault(:,1))+slipy.*cosd(fault(:,1)-90);
        
        slipx=slipnf(sizemat,slipx);
        slipy=slipnf(sizemat,slipy);
        slipe=slipnf(sizemat,slipe);
        slipn=slipnf(sizemat,slipn);
    end
    zslip(:,:,1) = slipx;
    zslip(:,:,2) = -slipy;
    return
else

slip=sum(reshape(X./nsub./1e6./3e10,[lensub,nsub]));
slip=slipnf(sizemat,slip);

stf=sum(sum(reshape(X,[lensub,sizemat]),2),3)*srate;

end

function slip=slipnf(sizemat,slip)

if sizemat==1
    slip=reshape(slip,sizemat);
    return;
else
    if size(sizemat,1)>1
        dex_sub=prod(sizemat');
        dex_sub=[0,cumsum(dex_sub)];
        slip0=slip;clear slip
        for i=1:size(sizemat,1)
            dex=dex_sub(i)+1:dex_sub(i+1);
            slip{i}=reshape(slip0(dex),sizemat(i,:));
        end
    end
end