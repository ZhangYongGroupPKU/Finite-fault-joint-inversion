function [slipg,Ggps,weit_gps,locagps] = getobg_gps(weight,obfile,sumob,...
    locasub,locadep,fault,grid,gridsize,trup_con)

fault = repmat(fault,[prod(grid),1]);

if weight
    load(obfile);
    locagps = gps(:,1:2);
    slipg = gps(:,3:5);
    weit_gps = weight*sqrt(sumob)/sqrt(slipg(:)'*slipg(:));
    [Ggps,~] = jt_gpsG(locagps,locasub,locadep,fault,gridsize,trup_con);
    Ggps = Ggps/3e16/prod(gridsize);
else
    slipg = [];
    Ggps = [];
    weit_gps = 0;
    locagps = [];
end
end