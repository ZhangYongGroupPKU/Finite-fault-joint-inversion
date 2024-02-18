function [G,sGDT,hdata,meangreen]=jt_seisG_moregreen(fault,trup_con0,delt, ...
    lensub,lambda1,lambda2,lambda3,obs,green,inter,srate,rake,sizemat)
%==========================================================================
% [X,rsq,G]=tempspac_invt0(fault,sizemat,gridsize,source,epi,...
%   lensub,lenrup,vrup,lambda1,lambda2,lambda3,obs,green,loca,phase,srate);
%  tempspac_invt is to inverse the spatio-temporal rupture process of
%  earthquake with multi-parameters.
%--------------------------------------------------------------------------
%Input
%       fault: each subfault has a fault mechanism,size: [nsub,3]
%        trup: in seconds
%    trup_con:
%        delt: unite is second

%Output:

%--------------------------------------------------------------------------
%                                                 Zhang Yong
%                                        2009/12/16 16:01   IGPCEA
%==========================================================================
%lensub0=lensub;
lensub=round(lensub*srate/inter);
%lenrup=round(lenrup./inter); % lenrup is unneccesary here because it has been
% used in the generation of trup_con0, and trup_con0 has been given for input parameters

meangreen=(max(abs(green2g(green,mean(fault)))));
%meangreen=(max(abs(green2g(green(:,:,:,1:42),mean(fault)))));
meangreen=mean(meangreen(:));
green=green/meangreen;
obs=obs/meangreen;

substf=tri_stf(inter,inter);substf(1)=[];
sg=size(green);
for i=1:sg(2)
    for j=1:sg(3)
        if length(sg)==3
            g0=conv(green(:,i,j),substf);
            green(:,i,j)=g0(inter:sg(1)+inter-1);
        else
            for k=1:sg(4)
                g0=conv(green(:,i,j,k),substf);
                green(:,i,j,k)=g0(inter:sg(1)+inter-1);
            end
        end
    end
end

delt=round(delt*srate);
%nsub=size(delt,1);
%srate1=srate/inter;
%[trup,trup_con0]=rup_get_trup(sizemat,source,gridsize,lensub,lenrup,vrup,srate1,tozero);

dex_out=stf2sub(1:length(trup_con0),size(delt,1),lensub);
% size(dex_out)
% size(trup_con0)
trup_con(dex_out)=trup_con0; %!


%trup_con=trup_con0(dex_out);
trup_con=trup_con(:);
xtrup_con=find(trup_con==0);
%ytrup_con=find(trup_con~=0);

% prepare for green's matrix:

G=conG_wch(green,delt,fault,lensub,trup_con,inter,rake);
%G=conG_wch(green,delt,fault,grid,ndot,trup_con,inter,rake)

% constrain in space:
if lambda1~=0
    
    D=sm_space_inter(sizemat,lensub);
    %D(xtrup_con(end:-1:1),:)=[];
    D(:,xtrup_con(end:-1:1))=[];

    sd=size(D);
    if rake==1
        D=[D,sparse(sd(1),sd(2));sparse(sd(1),sd(2)),D];
    end
end

% constrain in time:
if lambda2~=0
    [T]=sm_time_inter(sizemat,lensub);
    T(:,xtrup_con(end:-1:1))=[];
    %T(xtrup_con(end:-1:1),:)=[];
    
    sd=size(T);
    if rake==1
        T=[T,sparse(sd(1),sd(2));sparse(sd(1),sd(2)),T];
    end
end

% jointed matrix:
sGDT=[];
if lambda1~=0
    sGDT=[sGDT;D*lambda1];
    %clear D
end
if lambda2~=0
    sGDT=[sGDT;T*lambda2];
    %clear T
end
if lambda3~=0
    if rake~=1
        sGDT=sparse([sGDT;ones(1,sd(2))*lambda3]);
    else
        zlen=sd(2);
        sGDT=sparse([sGDT;[ones(1,zlen),zeros(1,zlen);zeros(1,zlen),ones(1,zlen)]*lambda3]);
    end
end


% data: obs and some zeros:

ssGDT=size(sGDT);
%hdata=[obs(:);zeros(ssGDT(1),1)];change
hdata=obs(:);
return
% solute the equation:
%[X,rsq]=cgls_xu_invt(full(GDT),hdata,1e-170,500,trup_con);

[X,rsq]=cgls_xu_invt0(G,sparse(sGDT),hdata,1e-170,500);
%[X,rsq]=cgls_xu_invt0(G,sparse(sGDT),hdata,1e-170,500);
Xinvt=X;

%rsq1=rsq*meangreen^2;
% syn0=G*X;
% rsq2=1-vr(obs(:),syn0(:));
%rsq=[min(rsq1),rsq2];

%[X,rsq]=cgls_zh0_spar(G,sparse(sGDT),hdata,100,'X(X<0)=0;');

%mean(abs(G(:)))
%XX=X;

syn=reshape(G*X*meangreen,size(obs));

G=G*meangreen;
%[X,rsq]=cgls_zh0_spar(full(GDT),sparse(sGDT),hdata,500,'X(X<0)=0;');
% X=fnnls(full([GDT;sGDT]),hdata);rsq=[];
% %==========-----------------------
if rake==1
    X0=zeros(lensub*nsub,2);
    X0(ytrup_con,:)=reshape(X,length(X)/2,2);
    X=X0;
else
    X0=zeros(lensub*nsub,1);
    X0(ytrup_con)=X;
    X=X0;
end
% %==========-----------------------

sX=size(X);
dex_out=sub2stf((1:sX(1))',nsub,round(lensub0*srate/inter));
if rake==0;
    XX(dex_out)=X;
else
    XX=zeros(length(dex_out),2);
    XX(dex_out,:)=X(:,:);
end
if inter>1
    if  rake==0;
        XX=reshape(XX,[round(lensub0*srate/inter),nsub]);
        XX=interm(XX,[inter,1]);
        sxx=size(XX);
        XX=[XX;zeros(lensub0*srate-sxx(1),sxx(2))];
        XX=XX(:);
    else
        XX1=reshape(XX(:,1),[round(lensub0*srate/inter),nsub]);
        XX1=interm(XX1,[inter,1]);
        XX2=reshape(XX(:,2),[round(lensub0*srate/inter),nsub]);
        XX2=interm(XX2,[inter,1]);
        
        sxx=size(XX1);
        XX1=[XX1;zeros(lensub0*srate-sxx(1),sxx(2))];
        XX2=[XX2;zeros(lensub0*srate-sxx(1),sxx(2))];
        XX=[XX1(:),XX2(:)];
    end
end
X=XX;
return
% ================================end======================================