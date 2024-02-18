function [X,rsq,qual]=cgls_xu_invt0(A,a,b,crit,itern)

%-----------------------------------------------------------------

% function [qual,X,rsq]=cgls_xu(A,X,b,Ac,g,h,crit2,bsq,...
%          rp,sdim,cont,velrupm,smthnss,itern);

% SPAR is a function used by SPARMAIN.

% For help see SPARMAIN.

%

% XL, SF 9-95

% Revised, Xu Lisheng, Paris, Aug., 1998
%-----------------------------------------------------------------                                                                                                                                                                                                                                                   

%residu=[];
%---------initial
sa=size(A);
Ac=sa(2);
crit2=Ac.*crit.^2;
X0=zeros(Ac,itern);
X=zeros(Ac,1);
bsq=sum(b.^2);
if isempty(a)
    xi = A*X-b;
else
    xi=[A*X;a*X]-b;
end
rp=sum(xi.^2);
if isempty(a)
    g=-A'*xi(1:sa(1));
else
    g=-A'*xi(1:sa(1))-a'*xi(sa(1)+1:end);
end
h=g;
residu=zeros(itern,1);
%-------------------------------
for iter=1:(10*Ac)
    if isempty(a)
        xi = A*h;
    else
    xi=[A*h;a*h];
    end
    anum=sum(g.*h);
    aden=sum(xi.^2);
    if aden==0
        disp('very singular matrix')
    end
    anum=anum./aden;
    %xi=X;
    X=X+anum.*h;
    %constrains: positive
    %constraints----------------------------------------------------------------------------------
%     xabs=abs(X);
%     X=(X+xabs)./2;

X(X<0)=0;

% X1=X(1:end-10);X2=X(end-9:end);
% X1(X1<0)=0;
% X=[X1;X2];
    
    %X(X<max(X)*0.01)=0;
    
    X0(:,iter)=X;
    %-----------------
    %X=X.*trup_con;


%     X=reshape(X,[lensub,grid_f(1),grid_f(2)]);
    %X(1,:,:)=0;
%     X(end,:,:)=0;
%     for i=1:grid_f(1)
%         for j=1:grid_f(2)
%             X(1:min(trup(i,j),end),i,j)=0;
%         end
%     end
%     X=X(:);
    % X(find(X<0))=0;
    %--------------------------------------------------------------------------------------------
    if isempty(a)
        xj = A*X - b;
    else
    xj=[A*X;a*X]-b;
    end
    rsq=sum(xj.^2);
    residu(iter)=rsq;
    if (rsq==rp||rsq<=(bsq*crit2)||iter==itern)
        %=================
        %here ITERN is the iteration time
        %qual=1;
        rsq=residu(1:iter);
        [xxx,numth]=min(rsq);
        X=X0(:,numth);
        %          disp('Solution is ready !')
        %          disp(['loop counter: ' num2str(iter)])
        return
    end
%     if rsq>rp 
%         X=xi;
%         qual=2;
%         disp('Wrong way,trying again')
%         %           disp(['for loop counter: ' num2str(iter)])
%         return
%     end
    %----------------
    rp=rsq;
    if isempty(a)
        xi=A'*xj(1:sa(1));
    else
    xi=A'*xj(1:sa(1))+a'*xj(sa(1)+1:end);
    end
    gg=sum(g.^2);
    dgg=sum((xi+g).*xi);
    %dgg=sum((xi).*xi);
    if gg==0
        qual=3;
        disp('gg=0')
        return
    end
    g=-xi;
    h=g+dgg./gg.*h;
end
qual=4;
rsq=residu;
disp('too many iterations')
%--------------------------End----------------------------------------