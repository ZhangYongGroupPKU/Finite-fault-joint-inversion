function dex_out=stf2sub(dex,sizemat,lensub)
%

dex_out=con2con(dex,prod(sizemat),lensub);


function mm0=con2con(mm,nsub,lensub)
%
sub_dex=ceil(mm./lensub);
dian_dex=mod(mm,lensub);

dian_dex(dian_dex==0)=lensub;

mm0=nsub*(dian_dex-1)+sub_dex;