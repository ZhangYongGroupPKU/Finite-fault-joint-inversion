function [slipg,Ggps,weit_gps,locagps] = getobg_gps2(weight,obfile,sumob,...
    locasub,locadep,fault,grid,gridsize,trup_con,rakevar)

fault = repmat(fault,[prod(grid),1]);

if weight
    load(obfile);
    locagps = gps(:,1:2);
    slipg = gps(:,3:5);
    weit_gps = weight*sqrt(sumob)/sqrt(slipg(:)'*slipg(:));
    if rakevar==1
        [Ggps,~] = jt_gpsG(locagps,locasub,locadep,fault,gridsize,trup_con);
    elseif rakevar==0
        [Ggps,~] = jt_gpsG2(locagps,locasub,locadep,fault,gridsize,trup_con);
    end
    Ggps = Ggps/3e16/prod(gridsize);
else
    slipg = [];
    Ggps = [];
    weit_gps = 0;
    locagps = [];
end
end