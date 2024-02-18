function [obs,dt,mco]=rup_obssyn_sm(obs,syn,srate)
%
%
if nargin==2
    srate=1;
end
obs = reshape(obs,length(obs),size(obs,2)/3,3);
syn = reshape(syn,length(syn),size(obs,2),3);
so=size(obs);

dt=zeros(so(2),1);
mco=zeros(so(2),1);

for i=1:so(2)
    [dt(i),mco(i)]=bestco_zh_sm(obs(:,i,:),syn(:,i,:));
end

xz=(abs(dt)>5*srate);
dt(xz)=0;mco(xz)=0;

for i=1:length(dt)
    if dt(i)>0
        obs(:,i,:)=[obs(dt(i)+1:end,i,:);zeros(dt(i),1,3)];
    elseif dt(i)<0
        obs(:,i,:)=[zeros(-dt(i),1,3);obs(1:end+dt(i),i,:)];
    end
end
obs = obs(:,:);