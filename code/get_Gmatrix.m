function [G,sGDT,hdata,meangreen,sumob,trup_con,locasub,locadep,trup] = get_Gmatrix(ob,g,...
    loca,epi,depth,fault,grid,gridsize,source,lensub,lenrup,vrup,srate,lam1,lam2,lam3,rakevar,set_edge0)

fault = repmat(fault,[prod(grid),1]);
lensub = lensub*srate;
lenrup = lenrup*srate;

if isempty(ob)
    loca = epi+[1,1];
end
[delt,locasub,locadep] = ruptime_nfnew(fault(1,:),grid([2,1]),...
    gridsize([2,1]),[1,source([2,1])],loca,epi,'P');
delt = zeros(size(delt));
locadep=locadep+depth-locadep(source(1));
[trup,trup_con0] = get_trup(grid,source,gridsize,lensub,lenrup,vrup,srate,set_edge0);
dex_out = stf2sub(1:length(trup_con0),grid,lensub);clear trup_con
trup_con(dex_out) = trup_con0; 
trup_con = trup_con(:);

if isempty(ob)
    G = [];sGDT = [];hdata = [];meangreen = 1;
    sumob = 1;
else
    [G,sGDT,hdata,meangreen] = jt_seisG_moregreen(fault,trup_con0(:),delt,...
    lensub/srate,lam1,lam2,lam3,ob,g,1,srate,rakevar,grid);
    sumob = hdata(:)'*hdata(:);
end
end
