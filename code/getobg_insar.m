function [slipsar,Gsar,weit_sar,locasar,nsar,dr] = getobg_insar(weight,obfile,sumob,...
    locasub,locadep,fault,grid,gridsize,trup_con)
fault = repmat(fault,[prod(grid),1]);
dr = [];
if weight
    load(obfile);
    locasar = insar(:,1:2);
    slipsar = insar(:,3);
    weit_sar = weight*sqrt(sumob)/sqrt(slipsar(:)'*slipsar(:));
    Gsar = jt_sarG(insar,locasub,locadep,prod(grid),locasar,fault,gridsize,trup_con);
    Gsar = Gsar/3e16/prod(gridsize);
else
    slipsar = [];
    Gsar = [];
    locasar = [];
    weit_sar = 0;
    nsar = [];
    dr = [];
end
end
