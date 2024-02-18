function G=geojt_conG(slipg,trup_con,israke)
%===============================================================================
%  G=geojt_conG(slipg,trup_con,israke);
%  Here it is a function to construct the matrix of green's function
%-------------------------------------------------------------------------------
%  Input:
%     slipg:  the green's function. has the size of [nsta,nsub] while israke==0. 
%             has the size of [nsta,nsub*2] while israke==1
%  trup_con:   
%    israke:  1 for rake variation, or else rake was fixed
%
%  Output:
%        G:  the matrix with the size of [nsta,length(find(trup_con==1))], or 
%            [nsta,length(find(trup_con==1))*2] while israke==1
%-------------------------------------------------------------------------------
%                                   Zhang Yong
%                                   2009/06/03  14:24     CEAIGP
%===============================================================================
if nargin==2
    israke=0;
end

ss=size(slipg);
nsub=ss(2);nsta=ss(1);
if israke==1
    nsub=nsub/2;
end
%lensub=length(trup_con)/nsub;

if israke==1
    %G=zeros(nsta,lensub*nsub*2);
    %G=zeros(nsta,length(find(trup_con==1))*2);
    trup_con=[trup_con;trup_con];
else
    %G=zeros(nsta,lensub*nsub);
    %G=zeros(nsta,length(find(trup_con==1)));
end
%sG=size(G);

% for i=1:nsta
%     submat=zeros(1,length(trup_con));
%     for j=1:nsub
%         if israke==1
%             submat(j:nsub:end/2)=slipg(i,j);
%             submat(j+length(submat)/2:nsub:end)=slipg(i,j+nsub);
%         else
%             submat(j:nsub:end)=slipg(i,j);
%         end
%     end
%     submat(xx)=[];
%     G(i,:)=submat;
% end

%-----------------

xx=trup_con==0;
onevec=ones(1,length(trup_con)/nsub/2);
submat=zeros(nsta,length(trup_con));
for j=1:nsub
    if israke==1
        submat(:,j:nsub:end/2)=slipg(:,j)*onevec;
%         submat(:,j+length(submat)/2:nsub:end)=slipg(:,j+nsub)*onevec;
        submat(:,j+size(submat,2)/2:nsub:end)=slipg(:,j+nsub)*onevec;
    else
        submat(:,j:nsub:end)=slipg(:,j)*ones(1,length(trup_con)/nsub);
    end
end
submat(:,xx)=[];
G=submat;
%=============================================end===================================================