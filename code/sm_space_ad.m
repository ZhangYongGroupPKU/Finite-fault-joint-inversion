function [x,mm,nn,num]=sm_space_ad(sizemat,lensub,flag)
%   clear
%   sizemat=[10,20];lensub=40;flag=1;
% =========================================================================
%  x=con_sm_mat_invt_multi0(sizemat,lensub,flag)
%
%  For smoothing the sliprate on the fault, x can be used in waveform
%  inversion of the rupture process, where lambda*x is a part of the
%  green's matrix, and lambda is a smoothing factor---whose value can be
%  determined by seeking the minimum of ABIC (Akaike's Bayesian Information
%  Criterion). (Akaike,1980)
%
%  sizemat is the nums of the fault's grids.sizemat(1) is for dip direction
%  and sizemat(2) for strike direction when dim==2. while dim==3,sizemat(1)
%  is for vertical direction, sizemat(2) for dip direction, and sizemat(3)
%  for strike direction.
%
%  lensub is the length of STFs of the sub-faults.
%
%  This program also can be used in moment tensor inversion study, where
%  lensub is the number of free parameters of the moment tensor for each
%  subfault.Usually lensub=6;
%
%  flag==1, edge elements donot change
%     else, edge elements are constrained to be zeros
%
%  See also con_sm_mat_invt, con_sm_mat, sm_mat, sm_mat1,
%  con_sm_mat_invt_multi.
%--------------------------------------------------------------------------
%                          Zhang Yong , Chen Yuntai and Xu Lisheng
%                             2006/09/25/   µØÕð¾ÖµØÇòËù (dim==2)
%                             2006/10/24/   Peking University (dim==3)
%==========================================================================
sb=sizemat;
dim=length(sb);
if nargin==2;flag=1;end
%x=sparse(prod(sb)*lensub,prod(sb)*lensub);

if dim==2
    mm=zeros(prod(sizemat)*5*lensub,1);
    nn=zeros(prod(sizemat)*5*lensub,1);
    num=zeros(prod(sizemat)*5*lensub,1);
    iter=0;
    
    nx=zeros(5,1);
    for i=1:sb(1)
        for j=1:sb(2)
            %-------------------------------------
            %          Strike: --------->        |
            %        Dip       2                 |
            %         |    1   3   5             |
            %         |        4                 |
            %        \|/                         |
            %-------------------------------------
            nx(1)=lensub*( (j-2)*sb(1)+ i -1 );
            nx(2)=lensub*( (j-1)*sb(1)+ i -1-1 );
            nx(3)=lensub*( (j-1)*sb(1)+ i -1 );
            nx(4)=lensub*( (j-1)*sb(1)+ i -1+1 );
            nx(5)=lensub*( (j  )*sb(1)+ i -1 );
            % find the positions of the four elements
            if flag==1
                %-------------------------------------------------
                % 4 edges:
                if i==1
                    nx(2)=NaN;
                end
                if i==sb(1)
                    nx(4)=NaN;
                end
                if j==1
                    nx(1)=NaN;
                end
                if j==sb(2)
                    nx(5)=NaN;
                end
                %--------------------------------------------------
                xn=find(nx>-1e0);
                for k=1:lensub
                    for zz=1:length(xn)
                        iter=iter+1;
                        if xn(zz)==3
                            mm(iter)=nx(3)+k;
                            nn(iter)=nx(3)+k;
                            num(iter)=-1;
                        else
                            mm(iter)=nx(3)+k;
                            nn(iter)=nx(xn(zz))+k;
                            num(iter)=1./(length(xn)-1);
                        end
                    end
                    %                     x(nx(3)+k,nx(xn)+k)=1./(length(xn)-1);
                    %                     x(nx(3)+k,nx(3)+k)=-1;
                end
            else
%                 for k=1:lensub
%                     % edges elements are set to be zeros
%                     if i==1||i==sb(1)||j==1||j==sb(2)
%                         x(nx(3)+k,nx(3)+k)=1e3;
%                     else
%                         x(nx(3)+k,nx([1,2,4,5])+k)=0.25;
%                         x(nx(3)+k,nx(3)+k)=-1;
%                     end
%                 end
            end % if or not : set edge elements to be zeros
        end % for i=1:sb(2)
    end % for i=1:sb(1)
elseif dim==3
    mm=zeros(prod(sizemat)*7*lensub,1);
    nn=zeros(prod(sizemat)*7*lensub,1);
    num=zeros(prod(sizemat)*7*lensub,1);
    iter=0;
    
    nx=zeros(7,1);
    for i=1:sb(1)
        for j=1:sb(2)
            for k=1:sb(3)
                % ---------------------------------------------------------
                %                           2      3
                %                           |    /
                %                           | /
                %                  1 ------ 4 ------- 7
                %                         / |
                %                      /    |
                %                    5      6
                %
                %                ___________________
                %               /                  /|
                %           i  /                  / |
                %             /__________________/  |
                %            |       k           |  |
                %          j |                   |  /
                %            |                   | /
                %            |___________________|/Left: the 3d fault model
                %
                %     3,5 ---- i th direction , the vertical direction
                %     2,6 ---- j th direction , also the dip direction
                %     1,7 ---- k th direction , also the strike direction
                %     nx(1) < nx(2) < nx(3) < nx(4) < nx(5) < nx(6) < nx(7)
                % ---------------------------------------------------------
                nx(1)=lensub*( (k-2)*sb(1)*sb(2) +(j-1)*sb(1) + i-1);
                nx(2)=lensub*( (k-1)*sb(1)*sb(2) +(j-2)*sb(1) + i-1);
                nx(3)=lensub*( (k-1)*sb(1)*sb(2) +(j-1)*sb(1) + i-1-1);
                nx(4)=lensub*( (k-1)*sb(1)*sb(2) +(j-1)*sb(1) + i-1);
                nx(5)=lensub*( (k-1)*sb(1)*sb(2) +(j-1)*sb(1) + i-1+1);
                nx(6)=lensub*( (k-1)*sb(1)*sb(2) +(j  )*sb(1) + i-1);
                nx(7)=lensub*( (k  )*sb(1)*sb(2) +(j-1)*sb(1) + i-1);
                if flag==1
                    %---------------------------------------------------
                    % smooth edge elements:
                    % 6 sides:
                    if i==1
                        nx(3)=NaN;
                    end
                    if i==sb(1)
                        nx(5)=NaN;
                    end
                    if j==1
                        nx(2)=NaN;
                    end
                    if j==sb(2)
                        nx(6)=NaN;
                    end
                    if k==1
                        nx(1)=NaN;
                    end
                    if k==sb(3)
                        nx(7)=NaN;
                    end
                    %-------------------------------------------------
                    xn=find(nx>-1e5);
                    for l=1:lensub
                        %                         x(nx(4)+l,nx(xn)+l)=1./(length(xn)-1);
                        %                         x(nx(4)+l,nx(4)+l)=-1;
                        for zz=1:length(xn)
                            if xn(zz)==4
                                iter=iter+1;
                                mm(iter)=nx(4)+l;
                                nn(iter)=nx(4)+l;
                                num(iter)=-1;
                            else
                                iter=iter+1;
                                mm(iter)=nx(4)+l;
                                nn(iter)=nx(xn(zz))+l;
                                num(iter)=1./(length(xn)-1);
                            end
                        end
                    end
                else
                    % set edges to be zeros
%                     for l=1:lensub
%                         if i==1||i==sb(1)||j==1||j==sb(2)||k==1||k==sb(3)
%                             % all sides: contain all edges and all corners
%                             x(nx(4)+l,nx(4)+l)=1e3;
%                         else
%                             x(nx(4)+l,nx(4)+l)=-1;
%                             x(nx(4)+l,nx([1,2,3,5,6,7])+l)=1/6;
%                         end
%                     end
                end % for flag
            end % for k=1:sb(3)
        end % for j=1:sb(2)
    end % for i=1:sb(1)
end % for dim

mm(iter+1:end)=[];
nn(iter+1:end)=[];
num(iter+1:end)=[];

x=sparse(mm,nn,num);
%================================end=======================================
%