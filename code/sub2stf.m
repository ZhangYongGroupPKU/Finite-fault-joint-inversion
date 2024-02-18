function dex_out=sub2stf(dex,sizemat,lensub)
% sizemat alos can be nsub
nsub=prod(sizemat);
dex_out=con2con(dex,nsub,lensub);

function mm0=con2con(mm,nsub,lensub)
%
sub_dex=ceil(mm./nsub);
dian_dex=mod(mm,nsub);

dian_dex(dian_dex==0)=nsub;

mm0=lensub*(dian_dex-1)+sub_dex;