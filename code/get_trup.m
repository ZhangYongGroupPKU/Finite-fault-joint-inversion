function [trup,trup_con]=get_trup(sizemat,source,gridsize,lensub,lenrup,vrup,srate,set_edge0)
%-------------------
trup=zeros(sizemat(1),sizemat(2));
trup_con=ones([lensub,sizemat]);
for i=1:sizemat(1) 
    for j=1:sizemat(2)
%         trup(i,j)=sqrt((abs((i-source(1))-0.5).*gridsize(1)).^2+...
%             (abs((j-source(2))-0.5).*gridsize(2)).^2)./vrup.*srate;
        trup(i,j)=sqrt((abs((i-source(1))).*gridsize(1)).^2+...
            (abs((j-source(2))).*gridsize(2)).^2)./vrup.*srate;
        trup_con(1:min(round(trup(i,j))+1,lensub),i,j)=0;
%         if i==sizemat(1)||i==1||j==1||j==sizemat(2)
%             trup_con(:,i,j)=0;%  edge to 0
%         end

        if set_edge0 == 0
            set0 = 0;
        elseif set_edge0 == 1
            set0 = (i==sizemat(1)||j==1||j==sizemat(2));
        elseif set_edge0 == 2
            set0 = (i==1||i==sizemat(1)||j==1||j==sizemat(2));
        else
            error('Wrong set_edge0!');
        end

        if set0
            trup_con(:,i,j)=0;%  edge to 0
        end

        trup_con(min(round(trup(i,j))+lenrup+2,lensub):end,i,j)=0;
    end 
end
trup_con=trup_con(:);