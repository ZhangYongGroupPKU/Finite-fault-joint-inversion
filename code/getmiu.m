function [miu1,miun] = getmiu(depth,fault,grid,gridsize,source,epi)

zeroloca = epi+[1,1];
[~,dep] = get_subloca(fault,grid([2,1]),gridsize([2,1]),[1,source([2,1])],zeroloca,epi);
dep=dep+depth-dep(source(1));%add
load('earth.mat');
miu = earth(:,4).*earth(:,3).^2.*1e9;
[~,y1] = histc(depth,earth(:,1));
[~,yn] = histc(dep(:),earth(:,1));
miu1 = miu(y1);
miun = miu(yn);
end
