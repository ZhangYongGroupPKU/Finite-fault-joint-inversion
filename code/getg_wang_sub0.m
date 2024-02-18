function [g_new]=getg_wang_sub0(fileg,leng,dist,t,t0,srate,dist0,rate)
%==========================================================================
%     [g_new]=getg_wang_sub(fileg,leng,dist,t,t0,srate,dist0)
%--------------------------------------------------------------------------
%  Input
%   fileg: the file of the green's function
%    leng: the length of the green's function
%    dist: epicentral distance, a vector
%       t: [t_p,t_s], the arrival times of the epicentral distances
%      t0: [t0_p,t0_s], the arrival times of the epicentral distances
%          contained in the database of green's functions
%   srate: sampling rate
%   dist0: the epicentral distances contained in the database of green's functions
%  Outpit 
%    g_new: obtained green's functions
%--------------------------------------------------------------------------
%    Zhang Yong, 2012-02-06 14:04, GFZ, Potsdam
%==========================================================================


% ensure the epicentral distance is monotonely inceased: it is not necessary
% because we 'sort' and 'unique' the index of the dist next
%[dist,dex]=sort(dist); 

% find the distance index
distnum=zeros(size(dist));
for i=1:length(dist)
    [~,nd]=min(abs(dist(i)-dist0));
    distnum(i)=nd;
end

% get the arrival times for the green's function in database
t0=t0(distnum,:);

% make unique and sort for the distance
[udnum,~,adex]=unique(distnum);


% the shift bytes for each epicentral distance
shift0=(leng+2)*10*4+4*3;
% the length of data read each time
leng0=(leng+2)*10+3;

ts=zeros(length(udnum),1);
g_temp=zeros(round(leng*rate/srate),10,length(udnum)); %3
%g_temp=zeros(leng,10,length(distnum)); %



fid=fopen(fileg,'r');
for i=1:length(udnum)
    % move the file pointer and read data
    fseek(fid,shift0*(udnum(i)-1),-1);
    z=fread(fid,leng0,'float');
        
    % the shift value
    ts(i)=z(2);
    %z(1:3)=[];
    z=z(4:end);
    
    % the data
    z=reshape(z,[leng+2,10]);
    %z([1,end],:)=[];
    z=z(2:end-1,:);
    
    % make resample for the green's functions and the sampling rate is
    % 'rate' now, not 'srate' again
    if rate~=srate
        z=resample(z,rate,srate);
    end
    
    g_temp(:,:,i)=z; % 1
%    g_temp(:,:,adex(i))=z; % 1+2+3
end
fclose(fid);

g_temp=g_temp(:,:,adex); %2
ts=ts(adex);


% shift the data by fixing their beginnings at the origin time
ts=round(ts*rate);
g=zeros(max(ts)+size(g_temp,1),size(g_temp,2),size(g_temp,3));
for i=1:size(g_temp,3)
    % get the synthetic seismograms by shift with ts
    if ts(i)<0
        zg=g_temp(1-ts(i):end,:,i);
        g(1:size(zg,1),:,i)=zg;
    else
%         zg=[zeros(ts(i),size(g_temp,2));g_temp(:,:,i)];
%         g(1:size(zg,1),:,i)=zg;
        g(ts(i)+1:ts(i)+size(g_temp,1),:,i)=g_temp(:,:,i);
    end
end

clear g_temp;

% make shifting for the P wave and S wave, respectively
t=round(t*rate);
t0=round(t0*rate);
tps=round(t(:,1)-t0(:,1));
tss=round(t(:,2)-t0(:,2));
g_new=zeros(size(g,1)+max(tss),size(g,2),size(g,3));

for i=1:size(g,3)
    gp=g(1:t(i,1),:,i);
    gps=g(t(i,1)+1:t(i,2),:,i);
    %[size(gp) tps(i)]
    gp_new=t_shift(gp,tps(i)); % for waves before P arrival times
    gps_new=t_shift(gps,tss(i)-tps(i)); % for waves between P ans S arrival times
    
    %g_new(1:size(g,1)+tss(i),:,i)=[gp_new;gps_new;g(t(i,2)+1:end,:,i)];
    g_new(1:size(gp_new,1),:,i)=gp_new;
    g_new(1+size(gp_new,1):size(gp_new,1)+size(gps_new,1),:,i)=gps_new;
    g_new(1+size(gp_new,1)+size(gps_new,1):size(g,1)+tss(i),:,i)=g(t(i,2)+1:end,:,i);
end
return
%-----------------------------------------------------------------------old


% while the dist is sorted
%g_new(:,:,dex)=g_new;