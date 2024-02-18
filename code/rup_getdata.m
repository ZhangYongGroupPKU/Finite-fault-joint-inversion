function [obs,mm,loca,chnl,finalpath,time]=rup_getdata(path,com,way,srate,d_range,epi)
%  [obs,mm,loca,epi,da,dep]=getdata_zh(path,com,way,srate)
%  get the SAC_ASC data from the files in 'datapath'
% Input:
%        path:   data folder
%        com:    such as 'Z' for vertical component, and 'A' for all components
%        way:    1 only get the data with the sampling rate of srate
%                0 get data with sampling rate more than srate   
%        srate:  sampling rate
% Output:
%        obs:  observed data
%        mm:   station codes
%        loca: location of stations
%        epi:  epicenter
%        da:   distance and azimuth
%

if isunix
    file=dir([path,'/*.SAC_ASC']);
else
    file1=dir([path,'\*.SAC_ASC']);
    file2=dir([path,'\*.SACA']);
    if length(file1)>length(file2)
        file=file1;
    else
        file=file2;
    end
end

nsta=length(file);
%prod([l,station_n])./1e8.*0.8


loca_ini=zeros(nsta,2);
num=0;
%------------------------------------------------------------------------

for n=1:nsta
    a=file(n).name;
    afile=[path,'/',a];
    filen=fopen(afile,'r');
    %headn=fscanf(filen,'%f');
    
    headn=zeros(110,1);
    for j=1:22
        heada=str2num(fgets(filen));
        headn((j-1)*5+1:j*5)=heada;
    end
    
    
        % way==1, only process the data with
    % a sampling rate equals srate:
    rate=round(1/headn(1));
    if rate<srate
        fclose(filen);
        continue
    end
    
    if headn(4)==-12345
        fclose(filen);
        continue
    end
    
    axx=repmat(' ',1e3,8);
    for i=1:8
        zaxx=fgets(filen);
        axx(1:length(zaxx),i)=zaxx;
    end

%     caxx=axx(1:2,2)';
%     if strcmp(caxx,'00')||strcmp(caxx,'  ')
%         %if caxx=='00'
%     else
%         fclose(filen);
%         continue
%     end
 
    %%判断分向，可以选择长周期或宽频带记录（第一个字符）。也可选择高增益或
    %%低增益记录（第二个字符）
    if strcmp(axx(17,7),'S');
        fclose(filen);
        continue
    end
    
%     caxxx=axx(18:18,7)';
%     if caxxx~='H'
%         fclose(filen);
%         continue
%     end                                
    
    %%去除掉L(B)H1和L(B)H2分量
    if axx(19,7)=='1'||axx(19,7)=='2'
        fclose(filen);
        continue
    end
    
    if strcmp(upper(com),'A')
    elseif ~strcmp(upper(com),axx(19,7))
        fclose(filen);
        continue
    end
    
    num=num+1;
    loca_ini(num,:)=[headn(32),headn(33)];
    usepath{num}=afile;


    
    fclose(filen);
end

if num==0
    error('Data sets wrong, can not read any data!!');
    return
end


loca_ini(num+1:end,:)=[];


da=distazim_zh(loca_ini,epi);
d_deg=da(:,1)./111.19492664455875;

% size(loca_ini)
% size(d_deg)
% size(usepath)
x1=find(d_deg<d_range(1));d_deg(x1)=[];loca_ini(x1,:)=[];usepath(x1')=[];
% size(loca_ini)
% size(d_deg)
% size(usepath)
x2=find(d_deg>d_range(2));d_deg(x2)=[];loca_ini(x2,:)=[];usepath(x2')=[];
% size(loca_ini)
% size(d_deg)
% size(usepath)

da=distazim_zh(loca_ini,epi);

% x3=azim_del(da(:,2),0.5);
% d_deg(x3)=[];loca_ini(x3,:)=[];usepath(x3')=[];

%------------------------------------------------------------------------


nsta=length(usepath);
if nsta==0
    error('Data sets wrong, can not read any data!!');
    return
end

leng=srate*5e3;

obs=zeros(leng,nsta);
mm=repmat(' ',nsta,5);
chnl=repmat(' ',nsta,6);
loca=zeros(nsta,2);
time=zeros(nsta,6);

num=0;
isepi=0;
for n=1:nsta
    filen=fopen(usepath{n},'r');
    %headn=fscanf(filen,'%f');
    
    headn=zeros(110,1);
    for j=1:22
        heada=str2num(fgets(filen));
        headn((j-1)*5+1:j*5)=heada;
    end
    
    % way==1, only process the data with
    % a sampling rate equals srate:
    rate=round(1/headn(1));
    if way==1
        if rate~=srate
            fclose(filen);
            continue
        end
        % otherwise, the sampling rate
        % should be more than required
    elseif rate<srate
        fclose(filen);
        continue
    end
    
    if headn(4)==-12345
        fclose(filen);
        continue
    end
    
    axx=repmat(' ',1e3,8);
    for i=1:8
        zaxx=fgets(filen);
        axx(1:length(zaxx),i)=zaxx;
    end

%     caxx=axx(1:2,2)';
%     if strcmp(caxx,'00')||strcmp(caxx,'  ')
%         %if caxx=='00'
%     else
%         fclose(filen);
%         continue
%     end
 
    %%判断分向，可以选择长周期或宽频带记录（第一个字符）。也可选择高增益或
    %%低增益记录（第二个字符）
%     caxxx=axx(18:18,7)';
%     if caxxx~='H'
%         fclose(filen);
%         continue
%     end
    
    %%去除掉L(B)H1和L(B)H2分量
    if axx(19,7)=='1'||axx(19,7)=='2'
        fclose(filen);
        continue
    end
    
    if strcmp(upper(com),'A')
    elseif ~strcmp(upper(com),axx(19,7))
        fclose(filen);
        continue
    end

    ob=fscanf(filen,'%f');
    fclose(filen);
    ob=ob(:);
    
    % clear the direct component:
    if length(ob)>100*rate
        ob=ob-mean(ob(1:100*rate));
    else
        continue
    end
    
    ob=ob/abs(headn(4));
%     if max(abs(ob))>0.15
%         continue
%     end
    
    % clear the incline component:
    %???
    

    if way~=1&&rate>=srate
        ob=resample(ob,srate,rate);
    end
        
    if headn(8)>0
        ob=ob(round(headn(8)*srate)+1:length(ob));
%        ob=ob(round(headn(8)*1)+1:length(ob));

    else
        ob=[zeros(abs(round(headn(8)*srate)),1);ob];
%        ob=[zeros(abs(round(headn(8)*1)),1);ob];
    end

    if max(abs(ob(1:min(leng,end))))==0
        continue
    end
    
    num=num+1;
    if length(ob)>leng
        obs(:,num)=ob(1:leng);
    else
        obs(1:length(ob),num)=ob;
    end
    
    
    % locations of stations
    loca(num,:)=[headn(32),headn(33)];

    mm(num,:)=axx(1:5,1)';
    
    chnl(num,:)=[axx(17:19,7)','_',axx(1:2,2)'];
    % only for the first time
    if num==1
        if headn(36)>-12345
            % there exist epi location
            isepi=1;
            epi=[headn(36),headn(37)];
        end
    end  
    
    finalpath{num}=usepath{n};
    
    time(num,:)=headn(71:76);
end


obs(:,num+1:end)=[];
mm(num+1:end,:)=[];
loca(num+1:end,:)=[];
chnl(num+1:end,:)=[];
time(num+1:end,:)=[];
time(:,5)=time(:,5)+time(:,6)./1e3;


disp([num2str(num) ' stations were selected successfully!!']);

% % there exist epi location
% if isepi==1
%     da=distazim_zh(loca,epi);
%     if strcmp(upper(com),'A')
%         da1=distazim_zh(epi,loca);
%         da=[da,da1(:,2)-180];
%     end
%     %%画出台站和震源的位置
%     epi_sta(epi,loca,mm,100)
%     camzoom(1.1)
% end
return