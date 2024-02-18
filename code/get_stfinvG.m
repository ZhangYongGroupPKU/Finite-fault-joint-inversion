function Gastf=get_stfinvG(grid,trup_con,tdt,lensub)
%%tdt 走时差 trup 起始破裂时间差,单位都是点
nsta=size(tdt,2);
nsub=grid(1)*grid(2);
tsf=tdt;   
Gastf=zeros(lensub*nsta,lensub*grid(1)*grid(2));
for i=1:nsub
    for j=1:nsta
        tem=zeros(1,lensub-abs(tsf(i,j)))+1;
        gtemp=diag(tem,-tsf(i,j));
        Gastf((j-1)*lensub+1:j*lensub,(i-1)*lensub+1:i*lensub)=gtemp;
    end
end
S=sparse(Gastf);
[m0,n0,v]=find(S);
m=con2con(m0,nsub,lensub);
n=con2con(n0,nsub,lensub);
Gastf=sparse(m0,n,v);
xtrup_con=find(trup_con==0);
Gastf(:,nsub*lensub)=0;
Gastf(:,xtrup_con(end:-1:1))=[];
end
function mm0=con2con(mm,nsub,lensub)
%
sub_dex=ceil(mm/lensub);
dian_dex=mod(mm,lensub);%

dian_dex(dian_dex==0)=lensub;

mm0=nsub*(dian_dex-1)+sub_dex;
end