function [obt,gt] = length_unification(obt,gt,obsm)
if isempty(obt) || isempty(obsm)
    return
else
    obt(size(obsm,1),:) = 0;
    gt(size(obsm,1),:,:,:) = 0;
end
end
