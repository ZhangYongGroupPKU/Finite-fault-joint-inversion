function [ob,g,mob] = normalize_obg(ob,g,obt,obsm,weight_tele,weight_sm)
mobt = sqrt(sum(obt.*obt));
if mobt == 0
    mobt = [];
end
mobsm = sqrt(sum(obsm.*obsm));
if mobsm == 0
    mobsm = [];
end
mob = [mobt/weight_tele,mobsm/weight_sm];
for i = 1:size(ob,2)
    ob(:,i) = ob(:,i)/mob(i);
    g(:,:,:,i) = g(:,:,:,i)/mob(i);
end
end
