function g=green2g(green,Mfp)
%==========================================================================
% green has the size of [leng,6,nsta] or [leng,6,nsta,nsub]
% Mfp can be the fault [strike dip rake], or the 6 elements of moment tensor
% g has the size of [leng,nsta] or [leng,nsta,nsub]
%==========================================================================
if length(Mfp)==3
    M=fp2mt_zh(1,Mfp(1),Mfp(2),Mfp(3));
else
    M=Mfp;
end
M=M(:);
if length(M)~=6
    error('elemets of moment tensor must be 6!');
end

sg=size(green);

if length(sg)==2
    g=green*M;
elseif length(sg)==3
    g=zeros(sg(1),sg(3));
    for i=1:sg(3)
        %                 for k=1:6
        %                     g(:,i)=g(:,i)+green(:,k,i)*M(k);
        %                 end
        g(:,i)=green(:,:,i)*M;
    end
else
    g=zeros(sg(1),sg(3),sg(4));
    for i=1:sg(3)
        for j=1:sg(4)
            %             for k=1:6
            %                 g(:,i,j)=g(:,i,j)+green(:,k,i,j)*M(k);
            %             end
            g(:,i,j)=green(:,:,i,j)*M;
        end
    end
end