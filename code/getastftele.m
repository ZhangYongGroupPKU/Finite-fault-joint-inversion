function [astftele,syn]=getastftele(obtele,g,source,grid,fault,lensub,band,srate)
g0=g(:,:,grid(1)*(source(2)-1)+source(1),:);
gastf=green2g(g0,fault);
gastf=gastf(:,:);
%[astftele,bt,~,cc]=pld(obtele,gastf,[lensub/10:lensub/10:lensub],1000,1);
[astftele,rsq,stfall,syn0]=pldf(obtele,gastf,lensub*srate,5000);
astftele=astftele*srate;
[bb,aa] = butter(3,band(2)*2/srate,'low');
astftele= filtfilt(bb,aa,astftele);
astftele(astftele<0)=0;
astftele=astftele(1:lensub*srate,:);
for i=1:size(gastf,2)
    syn(:,i)=conv(astftele(:,i),gastf(:,i))/srate;
end
syn(size(obtele)+1:end,:)=[];
end