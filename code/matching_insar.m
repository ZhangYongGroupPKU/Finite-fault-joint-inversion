function matching_insar(dlon,dlat,slipsar,synsar,locasar,nsar,dr)

mndeform = [-max(max(slipsar),-min(slipsar)),max(max(slipsar),-min(slipsar))];
for i = 1:length(nsar)
    if i==1
        dex = 1:nsar(1);
    else
        dex = nsar(i-1)+1:nsar(i);
    end
    
    loca = locasar(dex,:);
    
    if isempty(dr)
        numdata = min(2000,length(loca));
        loca1 = repmat(loca(1:numdata,:),numdata,1);
        loca2 = repmat_zh(loca(1:numdata,:),numdata);
        da = da_zh(loca1,loca2);
        dist = da(:,1);
        dist(1:numdata+1:end) = inf;
        r = min(dist)/2/111.1949;
    else
        r = dr(dex);
    end
    
    plt = {slipsar(dex),synsar(dex),slipsar(dex)-synsar(dex)};
    ttl = {'Observed','Synthetic','Residuals'};
        
    lontick = (ceil(min(loca(:,2))/dlon):floor(max(loca(:,2))/dlon))*dlon;
    lattick = (ceil(min(loca(:,1))/dlat):floor(max(loca(:,1))/dlat))*dlat;
    
    [lonlabel,latlabel] = tick2label(lontick,lattick);
    
    figure
    for j = 1:3
        subplot(1,3,j)
        plt_insardefor(loca,plt{j},mndeform,r,j);
        title(ttl{j});
        set(gca,'xtick',lontick);set(gca,'xticklabel',lonlabel);
        set(gca,'ytick',lattick);set(gca,'yticklabel',latlabel);
        set(gca,'fontsize',15);
    end
end
end
