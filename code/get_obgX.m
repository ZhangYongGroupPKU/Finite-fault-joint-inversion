function [X,ob,syn] = get_obgX(Xsolve,ob,syn,mob,lensub,trup_con,nsub,srate,rakevar)
lensub = lensub*srate;
if rakevar==1
X0 = zeros(lensub*nsub,2);
ytrup_con = find(trup_con ~= 0);
X0(ytrup_con,:) = reshape(Xsolve,length(Xsolve)/2,2);
dex_out = sub2stf((1:size(X0,1))',nsub,lensub);
X = zeros(length(dex_out),2);
X(dex_out,:) = X0(:,:);
else
X0 = zeros(lensub*nsub,1);
ytrup_con = find(trup_con ~= 0);
X0(ytrup_con) = Xsolve;
dex_out = sub2stf((1:size(X0,1))',nsub,lensub);
X = zeros(length(dex_out),1);
X(dex_out) = X0(:);
end
for i = 1:size(ob,2)
    ob(:,i) = ob(:,i).*mob(i);
    syn(:,i) = syn(:,i).*mob(i);
end
end
