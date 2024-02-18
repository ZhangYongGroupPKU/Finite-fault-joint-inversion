function [x]=sm_time_inter(grid,ndot)
% clear all
% grid=[8,15];
% ndot=60;

%=======================================
sizemat=grid;lensub=ndot;
%---------------------------------------------------------
m2=zeros(prod(sizemat)*lensub*3,1);
n2=zeros(prod(sizemat)*lensub*3,1);
num2=zeros(prod(sizemat)*lensub*3,1);
n=0;
nsub=prod(sizemat);
for i=1:lensub
    for j=1:nsub
        nx=(i-1)*nsub+j;
        
        if i==1
            n=n+1;
            m2(n)=nx;n2(n)=nx;num2(n)=-0.5;
            n=n+1;
            m2(n)=nx;n2(n)=nx+nsub;num2(n)=0.5;
            
        elseif i==lensub
            n=n+1;
            m2(n)=nx;n2(n)=nx;num2(n)=-0.5;
            n=n+1;
            m2(n)=nx;n2(n)=nx-nsub;num2(n)=0.5;
        else
            n=n+1;
            m2(n)=nx;n2(n)=nx;num2(n)=-1;
            n=n+1;
            m2(n)=nx;n2(n)=nx-nsub;num2(n)=0.5;
            n=n+1;
            m2(n)=nx;n2(n)=nx+nsub;num2(n)=0.5;
        end
    end
end
m2(n+1:end)=[];n2(n+1:end)=[];num2(n+1:end)=[];
x=sparse(m2,n2,num2);
%---------------------------------------------------------
% %----------------------------
% [x0,mm,nn,num]=sm_time_ad(sizemat,lensub);
% %[mm,nn,num]=sm_time_dex(sizemat,lensub);
% 
% mm0=con2con(mm,sizemat,lensub);
% nn0=con2con(nn,sizemat,lensub);
% x=sparse(mm0,nn0,num);
% %----------------------------------------------
% 
% function mm0=con2con(mm,sizemat,lensub)
% %
% sub_dex=ceil(mm/lensub);
% dian_dex=mod(mm,lensub);
% 
% dian_dex(dian_dex==0)=lensub;
% 
% mm0=prod(sizemat)*(dian_dex-1)+sub_dex;
% 
% 
% function [mm,nn,num]=sm_time_dex(sizemat,lensub)
% % ===================================================================================================
% % x=con_sm_mat_invt_time(sizemat,lensub,flag)
% % 
% %  for smoothing the sub-STFs of the subfaults on the fault, x can be used in waveform
% %  inversion of the rupture process, where lambda*x is a part of the green's
% %  matrix, and lambda is a smoothing factor---whose value can be determined by
% %  seeking the minimum of ABIC (Akaike's Bayesian Information Criterion).(Akaike,1980)
% %  
% %  sizemat is the nums of the fault's grids. sizemat(1) in dip direction
% %  and sizemat(2) in strike direction.
% %  lensub is the length of STFs of the sub-faults.
% %
% %  See also con_sm_mat_invt_multi,con_sm_mat_invt, con_sm_mat, sm_mat, sm_mat1
% %---------------------------------------------------------------------------------------------------
% %                                           Zhang Yong , Chen Yuntai and Xu Lisheng
% %                                                      2006/09/25/   地震局地球所
% %====================================================================================================
% sb=sizemat;
% dim=length(sb);
% x=sparse(prod(sb)*lensub,prod(sb)*lensub);
% 
% mm=zeros(prod(sizemat)*lensub*3,1);
% nn=zeros(prod(sizemat)*lensub*3,1);
% num=zeros(prod(sizemat)*lensub*3,1);
% n=0;
% if dim==2
%     for i=1:sb(1)
%         for j=1:sb(2)
%             for k=1:lensub
%                 nx=((j-1)*sb(1)+i-1)*lensub+k;
%                 if k==1
%                     n=n+1;
%                     mm(n)=nx;
%                     nn(n)=nx;
%                     num(n)=-0.5;
%                     n=n+1;
%                     mm(n)=nx;
%                     nn(n)=nx+1;
%                     num(n)=0.5;                    
%                     %x(nx,[nx,nx+1])=[-0.5,0.5];
%                 elseif k==lensub
%                     n=n+1;
%                     mm(n)=nx;
%                     nn(n)=nx-1;
%                     num(n)=0.5;
%                     n=n+1;
%                     mm(n)=nx;
%                     nn(n)=nx;
%                     num(n)=-0.5;                         
%                     %x(nx,[nx-1,nx])=[0.5,-0.5];
%                 else
%                     n=n+1;
%                     mm(n)=nx;
%                     nn(n)=nx-1;
%                     num(n)=0.5;
%                     n=n+1;
%                     mm(n)=nx;
%                     nn(n)=nx;
%                     num(n)=-1;
%                     n=n+1;
%                     mm(n)=nx;
%                     nn(n)=nx+1;
%                     num(n)=0.5;                    
%                     %x(nx,[nx-1,nx,nx+1])=[0.5,-1,0.5];
%                 end
%             end
%         end
%     end
%     mm(n+1:end)=[];
%     nn(n+1:end)=[];
%     num(n+1:end)=[];
%     
%     %x=sparse(mm,nn,num);
%     return
% elseif dim==3
%     % Many works need to do. It is so complex!
% end
% %================================end=======================================